function SmolPhi = PCEmaker(idx,qal,bounds,nvar,INFO)

Q = size(qal,1); %number of samples
SmolPhi = ones(size(qal,1),size(idx,1));

%scaling bound
for po = 1:nvar
    if INFO(po) == 1
        N{po}.Na = 2*((qal(:,po)-repmat(bounds(po,1),Q,1))./(repmat(bounds(po,2),Q,1)-repmat(bounds(po,1),Q,1)))-1;
    else
        N{po}.Na = (qal(:,po)-repmat(bounds(po,1),Q,1))./(sqrt(2*repmat(bounds(po,2),Q,1).^2));
    end
end

for ii = 1:size(idx,1) % Loop over PCE index
    h1 = ones(Q,1);
    ids = idx(ii,:);
    for jo = 1:length(ids) % Loop over number of index 2
        if INFO(jo) == 1
            %Legendre
            h1 = h1.* polyval(legendre(ids(jo),-1,1),N{jo}.Na)/sqrt(1/(2*ids(jo)+1));
        elseif INFO(jo) == 2
            %Hermite
            h1 = h1.* ((1/((2)^(ids(jo)/2)))*hermite(ids(jo),N{jo}.Na)/sqrt(factorial(ids(jo))));
        end
    end
     SmolPhi(:,ii) = h1;
end