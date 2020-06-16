
% Détection Linéaire zero-forcing 
%

function  dec = zf(r,H,M)
   
    G = (H'*H)\H';      % G = pseudo-inverse de H
    y = G*r;
    dec = qamdemod(y, M);
    
end

