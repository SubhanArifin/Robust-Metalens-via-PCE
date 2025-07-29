%% TOPOLOGY OPTIMIZATION OF METALENS (PCE)
clc; clear

%% INITIALIZATION
load('deterministic.mat','ResDVs')
realizations = 36;
lambda = lhsnorm(35,2.5^2,realizations);
eps_r  = 3*ones(realizations,1);

%% PRE-ALLOCATION
count = 0;
Densities = cell(1,length(0:0.25:7));
Data_ar = zeros(length(0:0.25:7),3);

%% PARETO OPTIMAL TRACING (POT)
tic;
for k = 0:0.25:7
    count = count + 1;
    [DVs,muFOM,sigFOM,~,loop,detMU,detSIG] = RobustAutolen_PCE(k,DVs,lambda,eps_r,realizations);
    Data_ar(count,:) = [muFOM sigFOM loop];
    Densities{1,count} = DVs;
    %load('deterministic.mat','DVs');  %ACTIVATE FOR MULTIPLE-RESTART
    DVs = ResDVs;
    if count == 1
        Data_det = [detMU detSIG];
    end
    clc
end
toc;