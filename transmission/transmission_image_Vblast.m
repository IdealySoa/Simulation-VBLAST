% Transmission d'image avec la technique V-BLAST
%% Clear
clear
clc
close all

%% Paramètres
Nt = 10;             
Nr = 10;             


ModType = 4;            % modulation : 2=BPSK, 4=QPSK, 16=16QAM, 64=64QAM
SNRindB = 20;

%% Lecture et traitement de l'image
t0 = tic;
img = imread('idealy-copie.jpg');
t = tic;
N = numel(img); % nombre de pixels
image_bin = de2bi(reshape(img,N,1));    % Répartition des éléments en une seule colonne et conversion en binaire
L = length(image_bin(:));           % La longueur de notre séquence de bit
txMsgBits = reshape(image_bin, log2(ModType), []);  % Ajustement des bits pour la modulation
txMapped = qammod(txMsgBits, ModType);  % Modulation 
x = reshape(txMapped,Nt,L/Nt);  % Arrangement des symboles pour les transmettre sur Nt antennes

% Transmission

SNR = 10^(SNRindB/10);
sigma = sqrt(1/SNR);
r_mmse_osic = zeros(L,length(SNRindB));

H = complex(rand(Nr,Nt),rand(Nr,Nt));  
r = awgn(H*x,SNRindB,'measured');                 % signal reçu

d_transmission = toc(t);
%% Réception
for col_idx = 1:L/Nt
    rsic = r(:,col_idx);
    r_mmse_osic((col_idx-1)*Nt+1:col_idx*Nt) = osic_mmse(rsic,H,sigma,ModType);
end

%% Reconstitution de l'image reçu
rec = mod(r_mmse_osic,2);
img_bin_rec = reshape(rec,size(image_bin));
img_rec = bi2de(img_bin_rec,2);
img_fin = uint8(reshape(img_rec,size(img)));    
d_execution_prog = toc(t0);
%% Affichage du résultat
figure
subplot(1,2,1)
imshow(img)
title('Image original')
subplot(1,2,2)
imshow(img_fin)
text = strcat('SNR=',num2str(SNRindB),',','VBLAST-MMSE/OSIC, Nt = ',num2str(Nt));
title({text;strcat(num2str(ModType),'-QAM')});
text_x = strcat('Durée d execution du programme',num2str(d_execution_prog),' s');
xlabel({strcat('Durée de transmission :',num2str(d_transmission),'s');text_x})



