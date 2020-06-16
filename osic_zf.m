% D�tection avec annulation d'interf�rence (ZF/OSIC)

% Matrice G : Nt*Nr
% Vecteur w : Nr*1
% Vecteur r : Nr*1

%%
function  rec = osic_zf(rsic,H,M)
    
    rec = zeros(1,length(rsic)); 

    [~,Nt] = size(H);
    k = zeros(1,Nt);    % Ordre de d�tection 

    for i = 1:Nt
        if(i==1)    % Initialisation 
            G = pinv(H);
            [~,k0]=min(sum(abs(G).^2,2));   % D�termination du premier symbole � �tre d�tecter 
            r = rsic;
        end
        
        k(i) = k0;
        y = G(k(i),:)*r;    % Estimation du symbole suivant l'ordre de d�tection d�fini par k(i)
        rec(k(i)) = qamdemod(y, M);
        a = qammod(rec(k(i)), M);
        r = r - a*H(:,k(i));    % Soustraction de la contribution de symbole(k(i))

        H(:,k(i)) = 0; % Remplacement du colonne de H(:,k(i)) par de zero car on a d�j� annul� la contribution de l'�metteur i
        G = pinv(H); % Nouveau pseudo-inverse
        temp = sum(abs(G).^2,2);  
        temp(k(1:i)) = Inf;
        [~,k0] = min(temp); % Nouvel �l�ment � d�tecter 
    end
end


