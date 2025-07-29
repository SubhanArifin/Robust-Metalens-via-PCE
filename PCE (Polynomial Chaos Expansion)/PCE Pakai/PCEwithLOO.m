function [bestPHI, bestOrder,LOOCV_err, R_squared] = PCEwithLOO(PMAX,x_exp,y_exp,boundst)

nvar    = size(x_exp,2); %random variabels dimension
pth     = 1; %polynomial trunctuation hyperbolic [0-1]
INFO    = 2*ones(1,nvar); %Legendre 1, Hermite 2

PHI = cell(PMAX); Err = zeros(PMAX,1);IDH = cell(PMAX);
for i = 1:PMAX
    order = i;
    min_samp = 2*factorial(nvar+order)/(factorial(nvar)*factorial(order));
    if size(x_exp,1) <= min_samp
        fprintf('For order %.2f, the number of basis (%.2f) is higher than number of samples (%.2f)\n',order,size(x_exp,1),min_samp)
        break
    end

    IDH{i} = totaltrunc(order,nvar,pth); %Generate Polynomial Index %Total ordere
    poly = PCEmaker(IDH{i},x_exp,boundst,nvar,INFO); %Generate Vandermonde Matrix %Experimental matrix
    PHI{i} = poly;
    
    % Residual computation NEW
    term1 = poly';
    term2 = poly'*poly;
    Mat = linsolve(term2,term1);
    coef = Mat*y_exp;
    pred = poly*coef;

    % LOOCV computation
    Ypred = zeros(1,size(x_exp,1));
    for j = 1:size(x_exp,1)
        x_temp = x_exp; x_temp(j,:)=[];
        y_temp = y_exp; y_temp(j,:)=[];

        TempPHI = PCEmaker(IDH{i},x_temp,boundst,nvar,INFO); %Generate Vandermonde Matrix %Experimental matrix
        Tempterm1 = TempPHI'*y_temp;
        Tempterm2 = TempPHI'*TempPHI;
        Tempcoef = linsolve(Tempterm2,Tempterm1);

        TempPHI2 = PCEmaker(IDH{i},x_exp(j,:),boundst,nvar,INFO);
        Ypred(j) = TempPHI2*Tempcoef;
    end

    %% Residual
    resd = (pred - Ypred') ./ diag(poly*Mat);
    Err(i) = 1/size(y_exp,1) * sum(( resd.^2 ));
end
    [LOOCV_err,indx] = min(Err);
    bestOrder = indx;
    bestPHI = PHI{indx};
    R_squared = 1-LOOCV_err/var(y_exp); %determination coefficient
   
