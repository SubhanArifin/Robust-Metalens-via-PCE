function [ery, erloo, erlon, yloo] = LOOcalc(PHIS,ypr,ytr)
%PHIE = PHIS*inv(PHIS'*PHIS)*PHIS';
PHIE=PHIS*((PHIS'*PHIS)\PHIS');
dp = diag(PHIE);
% dp
% pause
erloo = (ytr-ypr)./(1-dp(1:length(ypr)));
%ery = sqrt(mean(erloo.^2));
ery = mean(abs(erloo));
%ery = (ery).^0.5; %kemas edit
msigma = mean(ypr);
erlon = mean(erloo.^2)./mean((ytr-msigma).^2);
yloo = ytr-erloo;
