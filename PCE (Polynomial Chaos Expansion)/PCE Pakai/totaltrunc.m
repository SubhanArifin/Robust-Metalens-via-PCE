function idx = totaltrunc(nix,d,q)
% inputs:
% nix = Maximum polynomial order
% d = number of variables
% q = hyperbolic truncation parameter

% Generate index for polynomial chaos expansion
idx = [];
for i = nix:-1:0
    idx = [idx;MonCof(i,d)];
end
idx = flipud(idx);

% Now truncate further! (if q == 1, this equals to total-order expansion)
if q<1
idp = sum(idx.^q,2).^(1/q);
idx = idx(idp <=(nix+0.000001),:);
end
