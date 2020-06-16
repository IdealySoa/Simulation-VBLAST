% Détection avec annulation d'interférence (MMSE/OSIC)



function  rec = osic_mmse(rsic,H,sigma,M)
    
    rec = zeros(1,length(rsic));

    [~,Nt] = size(H);
    k = zeros(1,Nt); % k est l'ordre de détection

    for i = 1:Nt
        % Initalisation
        if(i==1)
            G = (H'*H+sigma.^2*eye(Nt))\H'; 
            [~,k0]=min(sum(abs(G).^2,2));   % Détermination du premier symbole à détecter
            r = rsic;
        end
        
        k(i) = k0;
        y = G(k(i),:)*r;        % Estimation du symbole suivant l'ordre de détection défini par k(i)        
        rec(k(i)) = qamdemod(y, M);
        a = qammod(rec(k(i)), M);
        r = r- a*H(:,k(i));     % Annulation de l'effet du symbole k(i)

        H(:,k(i)) = 0;          % Réduction de la taille de la matrice canal
        G = (H'*H+sigma.^2*eye(Nt))\H';
        temp = sum(abs(G).^2,2);
        temp(k(1:i)) = Inf;
        [~,k0] = min(temp);     % Détermination du nouvel élement à etre détecté
    end
end