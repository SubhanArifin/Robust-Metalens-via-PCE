function PHI = PCE(PMAX,x_exp,boundst)

nvar    = size(x_exp,2); %random variabels dimension
pth     = 1; %polynomial trunctuation hyperbolic [0-1]
INFO    = 1*ones(1,nvar); %Legendre 1, Hermite 2

IDH = totaltrunc(PMAX,nvar,pth); %Generate Polynomial Index
PHI = PCEmaker(IDH,x_exp,boundst,nvar,INFO); %Generate Vandermonde Matrix %Experimental matrix