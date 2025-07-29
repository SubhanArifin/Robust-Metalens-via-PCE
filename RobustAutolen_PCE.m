function [ResDVs,muFOM,sigFOM,Objeval,loop,detFOM_mu,detFOM_sig] = RobustAutolen_PCE(k,DVs,lambda,eps_r,realizations)
%% INITIALIZATION
nElX= 400; nElY = 200; DesignThicknessElements = 15;
DDIdx = repmat([1:nElY:nElX*nElY],DesignThicknessElements,1);
dVElmIdx = DDIdx+repmat([165:165+DesignThicknessElements-1]',1,nElX);
targetXY=[200,80]; fR=6.0; MaxItr=200;

%% PROPERTIES
% SETUP OF PHYSICS PARAMETERS
phy.scale = 1e-9; % Scaling finite element side length to nanometers

% SETUP OF ALL INDEX SETS, ELEMENT MATRICES AND RELATED QUANTITIES
dis.nElX = nElX; % number of elements in x direction
dis.nElY = nElY; % number of elements in y direction
dis.tElmIdx = (targetXY(1)-1)*nElY+targetXY(2); % target index
dis.dVElmIdx = dVElmIdx; % design field element indices in model of physics
[dis.LEM,dis.MEM] = ELEMENT_MATRICES(phy.scale);
[dis]=INDEX_SETS_SPARSE(dis); % Index sets for discretized model

% SETUP FILTER AND THRESHOLDING PARAMETERS
filThr.beta = 5; % Thresholding sharpness
filThr.eta = 0.5; % Thresholding level
[filThr.filKer, filThr.filSca] = DENSITY_FILTER_SETUP( fR, nElX, nElY);

% INITIALIZE DESIGN VARIABLES, BOUNDS AND OPTIMIZER OPTIONS
%dVs_in = ones(length(dis.dVElmIdx(:)),1)*dVini; % Design variables
dVs_in = DVs;

% DISTRIBUTE MATERIAL IN MODEL DOMAIN BASED ON DESIGN FIELD
dFP(1:dis.nElY,1:dis.nElX) = 0; % Design field in physics, 0: air
dFP(dis.nElY:-1:ceil(dis.nElY*9/10),1:dis.nElX) = 1; % 1: material

%% INITIALIZE ITERATION 
LBdVs = zeros(length(dVs_in),1); % Lower bound on design variables
UBdVs = ones(length(dVs_in),1); % Upper bound on design variables
dVs     = [0.9*dVs_in 0.5*dVs_in dVs_in zeros(length(dis.dVElmIdx(:)),MaxItr)];
temp    = 0.5*(max(dVs_in)-min(dVs_in));
low     = dVs_in - temp; %initlow
upp     = dVs_in + temp; %initup

loop       = 0;
r_ch       = 1;
m_value    = zeros(1,MaxItr);
s_value    = zeros(1,MaxItr);
F_value    = zeros(1,MaxItr);
chg        = zeros(1,MaxItr);
mat        = cell(1,MaxItr);
Elf        = cell(1,MaxItr);

while r_ch > 1e-6
    loop = loop+1;
    %% MMA X VALUE
    dVs_val    = dVs(:,loop+2);
    xold1      = dVs(:,loop+1);
    xold2      = dVs(:,loop);
   
    dFP(dis.dVElmIdx(:)) = dVs_val; % Design variables inserted in design field
    
    % FILTERING THE DESIGN FIELD AND COMPUTE THE MATERIAL FIELD
    dFPS = DENSITY_FILTER(filThr.filKer,filThr.filSca,dFP,ones(dis.nElY,dis.nElX));
    dFPST = THRESHOLD( dFPS, filThr.beta, filThr.eta);

    %% RANDOM REALIZATIONS
    FOM_f = zeros(1,realizations);
    sensFOM_f= zeros(DesignThicknessElements,nElX,realizations);
    for     i = 1:realizations    
        % COMPUTE THE MATERIAL FIELD
        phy.eps_r = eps_r(i);    
        phy.k = 2*pi/(lambda(i)*phy.scale);
        [A,dAdx] = MATERIAL_INTERPOLATION(phy.eps_r,dFPST,1.0); % Material field
        
        % CONSTRUCT THE SYSTEM MATRIX
        [dis,F] = BOUNDARY_CONDITIONS_RHS(phy.k,dis,phy.scale);
        dis.vS = reshape(dis.LEM(:)-phy.k^2*dis.MEM(:)*(A(:).'),16*dis.nElX*dis.nElY,1);
        S = sparse([dis.iS(:);dis.iBC(:)],[dis.jS(:);dis.jBC(:)],[dis.vS(:);dis.vBC(:)]);
            
        % SOLVING THE STATE SYSTEM: S * Ez = F
        [L,U,Q1,Q2] = lu(S); % LU - factorization
        Ez = Q2 * (U\(L\(Q1 * F))); Ez = full(Ez); % Solving
            
        % FIGURE OF MERIT
        P = sparse(dis.edofMat(dis.tElmIdx,:),dis.edofMat(dis.tElmIdx,:),1/4,...
        (dis.nElX+1)*(dis.nElY+1),(dis.nElX+1)*(dis.nElY+1)); % Weighting matrix
        FOM = Ez' * P * Ez; % Solution in target element
        
        % ADJOINT RIGHT HAND SIDE
        AdjRHS = P*(2*real(Ez) - 1i*2*imag(Ez));
        
        % SOLVING THE ADJOING SYSTEM: S.’ * AdjLambda = AdjRHS
        AdjLambda = (Q1.') * ((L.')\((U.')\((Q2.') * (-1/2*AdjRHS)))); % Solving
        
        % COMPUTING SENSITIVITIES
        dis.vDS = reshape(-phy.k^2*dis.MEM(:)*(dAdx(:).'),16*dis.nElX*dis.nElY,1);
        DSdx = sparse(dis.iElFull(:),dis.jElFull(:),dis.vDS(:)); % Constructing dS/dx
        DSdxMulV = DSdx * Ez(dis.idxDSdx); % Computing dS/dx * Field values
        DsdxMulV = sparse(dis.iElSens,dis.jElSens,DSdxMulV);
        sens = 2*real(AdjLambda(dis.idxDSdx).' * DsdxMulV); % Computing sensitivites
        sens = full(reshape(sens,dis.nElY,dis.nElX));
        
        % FILTERING SENSITIVITIES
        DdFSTDFS = DERIVATIVE_OF_THRESHOLD( dFPS, filThr.beta, filThr.eta);
        sensFOM = DENSITY_FILTER(filThr.filKer,filThr.filSca,sens,DdFSTDFS);
        
        % EXTRACTING SENSITIVITIES FOR DESIGNABLE REGION
        sensFOM = sensFOM(dis.dVElmIdx); %15x400
        
        % FMINCON DOES MINIMIZATION
        FOM_f(i) = FOM; 
        sensFOM_f(:,:,i) = sensFOM;

        % PRINT REALIZATION RESULTS
        fprintf('k:%.2f  It:%5i  realizations:%6i\n',k,loop,i);
    end
    %% POLYNOMIAL CHAOS EXPANSION (PCE)
    PHI = PCE(5,[lambda],[35 2.5]); %Add PCE folder to path / working directory
    % OBJECTIVE FUNCTION 'F' (FOM)
    COFO          = PHI\FOM_f'; %Coef c using OLS % c = FOM_f
    m_value(loop) = COFO(1); %Expected value of compliance
    s_value(loop) = sqrt(sum(COFO(2:end).^2)); %Std value of compliance
    F_value(loop) = -(m_value(loop) - k*s_value(loop)); %Objective Function
    % RANDOM RESPONSE SENSITIVITY 'dF' ANALYSIS
    dCOFO         = zeros(DesignThicknessElements,nElX,length(COFO));
    for ii=1:DesignThicknessElements
        for jj=1:nElX
            dCOFO(ii,jj,:) = PHI\reshape(sensFOM_f(ii,jj,:),realizations,1);
        end
    end
    %Mean Sensitivity
    dm = dCOFO(:,:,1);
    %Std Sensitivity
    ds = zeros(DesignThicknessElements,nElX);
    for ii=2:size(COFO,1)
        ds = ds + COFO(ii)*dCOFO(:,:,ii)/s_value(loop);
    end
    %Objective Function Sensitivity
    dF = -(dm - (k*ds));
    %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % Design Update
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(1,length(dis.dVElmIdx(:)),loop,dVs_val,LBdVs,UBdVs,xold1,xold2,...
        F_value,dF(:),0,0,low,upp,1,0,1e3,0);
    dVs(:,loop+3) = xmma;
    change = max(abs(dVs_val(:) - xold1(:)));
    chg(loop) = change;
    if loop > 10
        r_ch = abs(chg(loop-1)-chg(loop));
    end
    %% PRINT RESULTS
    fprintf('It.:%d  Obj.:%.3f FOM:%.3f std.:%.3f R.ch.:%.6f\n',loop,-F_value(loop),m_value(loop),s_value(loop),r_ch);
    %% PLOT DENSITIES
    %figure(1); % Field intensity, |Ez|^2
    %imagesc((reshape(Ez.*conj(Ez),dis.nElY+1,dis.nElX+1))); colorbar; axis equal;
    %figure(2); % Physical design field
    %imagesc(1-dFPST); colormap(gray); axis equal; drawnow;
    %% STOPPING CRITERIA
    if loop > MaxItr - 1
        break
    end
end 
%% DATA COLLECTION
muFOM     = m_value(loop);
sigFOM    = s_value(loop);

detFOM_mu  = m_value(1);
detFOM_sig = s_value(1);

Objeval   = -F_value(loop);
ResDVs    = dVs(:,loop+3);
end 

%% AUXILIARY FUNCTIONS 
%%%%%%%%%%%% ABSORBING BOUNDARY CONDITIONS AND RIGHT HAND SIDE %%%%%%%%%%%%
function [dis,F] = BOUNDARY_CONDITIONS_RHS(waveVector,dis,scaling)
AbsBCMatEdgeValues = 1i*waveVector*scaling*[1/6 ; 1/6 ; 1/3 ; 1/3];
% ALL BOUNDARIES HAVE ABSORBING BOUNDARY CONDITIONS
dis.iBC = [dis.iB1(:);dis.iB2(:);dis.iB3(:);dis.iB4(:)];
dis.jBC = [dis.jB1(:);dis.jB2(:);dis.jB3(:);dis.jB4(:)];
dis.vBC = repmat(AbsBCMatEdgeValues,2*(dis.nElX+dis.nElY),1);
% BOTTOM BOUNDARY HAS INCIDENT PLANE WAVE
F = zeros((dis.nElX+1)*(dis.nElY+1),1); % System right hand side
F(dis.iRHS(1,:)) = F(dis.iRHS(1,:))+1i*waveVector;
F(dis.iRHS(2,:)) = F(dis.iRHS(2,:))+1i*waveVector;
F = scaling*F;
end
%%%%%%%%%%%%%%%%%%%%% CONNECTIVITY AND INDEX SETS %%%%%%%%%%%%%%%%%%%%%%%%%
function [dis]=INDEX_SETS_SPARSE(dis)
% INDEX SETS FOR SYSTEM MATRIX
nEX = dis.nElX; nEY = dis.nElY; % Extracting number of elements
nodenrs = reshape(1:(1+nEX)*(1+nEY),1+nEY,1+nEX); % Node numbering
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nEX*nEY,1); % First DOF in element
dis.edofMat = repmat(edofVec,1,4)+repmat([0 nEY+[1 0] -1],nEX*nEY,1);
dis.iS = reshape(kron(dis.edofMat,ones(4,1))',16*nEX*nEY,1);
dis.jS = reshape(kron(dis.edofMat,ones(1,4))',16*nEX*nEY,1);
dis.idxDSdx = reshape(dis.edofMat',1,4*nEX*nEY);
% INDEX SETS FOR BOUNDARY CONDITIONS
TMP = repmat([[1:nEY];[2:nEY+1]],2,1);
dis.iB1 = reshape(TMP,4*nEY,1); % Row indices
dis.jB1 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEY,1); % Column indices
TMP = repmat([1:(nEY+1):(nEY+1)*nEX;(nEY+1)+1:(nEY+1):(nEY+1)*nEX+1],2,1);
dis.iB2 = reshape(TMP,4*nEX,1);
dis.jB2 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEX,1);
TMP = repmat([(nEY+1)*(nEX)+1:(nEY+1)*(nEX+1)-1;(nEY+1)*(nEX)+2:(nEY+1)*(nEX+1)],2,1);
dis.iB3 = reshape(TMP,4*nEY,1);
dis.jB3 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEY,1);
TMP = repmat([2*(nEY+1):nEY+1:(nEY+1)*(nEX+1);(nEY+1):nEY+1:(nEY+1)*(nEX)],2,1);
dis.iB4 = reshape(TMP,4*nEX,1);
dis.jB4 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEX,1);
dis.iRHS = TMP;
% INDEX SETS FOR INTEGRATION OF ALL ELEMENTS
ima0 = repmat([1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4],1,nEX*nEY).';
jma0 = repmat([1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4],1,nEX*nEY).';
addTMP = repmat(4*[0:nEX*nEY-1],16,1);
addTMP = addTMP(:);
dis.iElFull = ima0+addTMP;
dis.jElFull = jma0+addTMP;
% INDEX SETS FOR SENSITIVITY COMPUTATIONS
dis.iElSens = [1:4*nEX*nEY]';
jElSens = repmat([1:nEX*nEY],4,1);
dis.jElSens = jElSens(:);
end
%%%%%%%%%%%%%%%%%% MATERIAL PARAMETER INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%
function [A,dAdx] = MATERIAL_INTERPOLATION(eps_r,x,alpha_i)
A = 1 + x*(eps_r-1) - 1i * alpha_i * x .* (1 - x); % Interpolation
dAdx = (eps_r-1)*(1+0*x) - 1i * alpha_i * (1 - 2*x); % Derivative of interpolation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xS]=DENSITY_FILTER(FilterKernel,FilterScaling,x,func)
xS = conv2((x .* func)./FilterScaling,FilterKernel,'same');
end
function [ Kernel, Scaling ] = DENSITY_FILTER_SETUP( fR, nElX, nElY )
[dy,dx] = meshgrid(-ceil(fR)+1:ceil(fR)-1,-ceil(fR)+1:ceil(fR)-1);
Kernel = max(0,fR-sqrt(dx.^2+dy.^2)); % Cone filter kernel
Scaling = conv2(ones(nElY,nElX),Kernel,'same'); % Filter scaling
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ xOut ] = THRESHOLD( xIn, beta, eta)
xOut = (tanh(beta*eta)+tanh(beta*(xIn-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));
end
function [ xOut ] = DERIVATIVE_OF_THRESHOLD( xIn, beta, eta)
xOut = (1-tanh(beta*(xIn-eta)).^2)*beta./(tanh(beta*eta)+tanh(beta*(1-eta)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ELEMENT MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LaplaceElementMatrix,MassElementMatrix] = ELEMENT_MATRICES(scaling)
% FIRST ORDER QUADRILATERAL FINITE ELEMENTS
aa=scaling/2; bb=scaling/2; % Element size scaling
k1=(aa^2+bb^2)/(aa*bb); k2=(aa^2-2*bb^2)/(aa*bb); k3=(bb^2-2*aa^2)/(aa*bb);
LaplaceElementMatrix = [k1/3 k2/6 -k1/6 k3/6 ; k2/6 k1/3 k3/6 -k1/6; ...
-k1/6 k3/6 k1/3 k2/6; k3/6 -k1/6 k2/6 k1/3];
MassElementMatrix = aa*bb*[4/9 2/9 1/9 2/9 ; 2/9 4/9 2/9 1/9 ; ...
1/9 2/9 4/9 2/9; 2/9 1/9 2/9 4/9];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% FOM EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FOM] = OBJECTIVE_eval(dVs,dis,phy,filThr)
% DISTRIBUTE MATERIAL IN MODEL DOMAIN BASED ON DESIGN FIELD
dFP(1:dis.nElY,1:dis.nElX) = 0; % Design field in physics, 0: air
dFP(dis.nElY:-1:ceil(dis.nElY*9/10),1:dis.nElX) = 1; % 1: material
dFP(dis.dVElmIdx(:)) = dVs; % Design variables inserted in design field

% FILTERING THE DESIGN FIELD AND COMPUTE THE MATERIAL FIELD
dFPS = DENSITY_FILTER(filThr.filKer,filThr.filSca,dFP,ones(dis.nElY,dis.nElX));
dFPST = THRESHOLD( dFPS, filThr.beta, filThr.eta);
[A,~] = MATERIAL_INTERPOLATION(phy.eps_r,dFPST,1.0); % Material field

% CONSTRUCT THE SYSTEM MATRIX
[dis,F] = BOUNDARY_CONDITIONS_RHS(phy.k,dis,phy.scale);
dis.vS = reshape(dis.LEM(:)-phy.k^2*dis.MEM(:)*(A(:).'),16*dis.nElX*dis.nElY,1);
S = sparse([dis.iS(:);dis.iBC(:)],[dis.jS(:);dis.jBC(:)],[dis.vS(:);dis.vBC(:)]);
%tic;

% SOLVING THE STATE SYSTEM: S * Ez = F
[L,U,Q1,Q2] = lu(S); % LU - factorization
Ez = Q2 * (U\(L\(Q1 * F))); Ez = full(Ez); % Solving
%toc;

% FIGURE OF MERIT
P = sparse(dis.edofMat(dis.tElmIdx,:),dis.edofMat(dis.tElmIdx,:),1/4,...
(dis.nElX+1)*(dis.nElY+1),(dis.nElX+1)*(dis.nElY+1)); % Weighting matrix
FOM = Ez' * P * Ez; % Solution in target element

% PLOTTING AND PRINTING
figure(1); % Field intensity, |Ez|^2
imagesc((reshape(Ez.*conj(Ez),dis.nElY+1,dis.nElX+1))); colorbar; axis equal;
figure(2); % Physical design field
imagesc(1-dFPST); colormap(gray); axis equal; drawnow;
disp(['FOM: ' num2str(FOM)]); % Display FOM value
end