% Détection linéaire MMSE 


function  dec = mmse(r,H,sigma,M)
    [~,Nt] = size(H);
    W_mmse = (H'*H+sigma.^2*eye(Nt))\H';
    y = W_mmse*r;
    dec = qamdemod(y, M);
    
end


