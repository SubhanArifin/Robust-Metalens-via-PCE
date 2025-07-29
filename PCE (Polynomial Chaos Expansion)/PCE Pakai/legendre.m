function La = legendre(n,lb,ub)

% check n
if( n<0 ), error('The order of legendre polynomial must be greater than or equal to 0.'); end

% again check n is an integer
if( 0~=n-fix(n) ), error('The order of legendre polynomial must be an integer.'); end

% call the legendre recursive function.
L0 = 1;
L1 = [1 -((ub+lb)/2)];

if n == 0
    La = L0;
    return
elseif n  == 1
    La = L1;
    if lb == 0
         La = [2 -1];
    end
    return
end
        

% Perform Gram Schmidt orthogonalization
for i = 1:n-1
    if i==1;K = L0;L=L1;end
    
    a = (polyval(polyint(conv([1 0],conv(L,L))),ub)-polyval(polyint(conv([1 0],conv(L,L))),lb)) / (polyval(polyint(conv(L,L)),ub)- polyval(polyint(conv(L,L)),lb));  
    h1 = conv(polymin([1 0],a),L);
    
    b = (polyval(polyint(conv(L,L)),ub)-polyval(polyint(conv(L,L)),lb))/(polyval(polyint(conv(K,K)),ub)-polyval(polyint(conv(K,K)),lb));
    h2 = conv(b,K);
    
    La = polymin(h1,h2);
    
    K = L;
    L = La;  
end

%Perform normalization if needed
cons = 1/polyval(La,1);
La = La*cons;


           
    




