function y = PCEpredND(qal,coff,idx,bounds,nvar,INFO)

Q = size(qal,1);
% Q
PHIP = PCEmaker(idx,qal,bounds,nvar,[INFO]);
 
% size(coff)
y = (PHIP*coff);
