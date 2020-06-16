%% Simulation VBLAST : simulation.m
clear;
clc;
close all;

%% Paramètres
Nt = 4;             % Nombre d'antennes à l'émission
Nr = 4;             % Nombre d'antennes à la reception

L = 100;            % Nombre de symboles par trame
Reps = 50;     % Boucle de répétition par SNR

SNRindB = 0:2:40;       % Valeur SNR
M = 64;             % modulation : 4=QPSK, 16=16QAM, 64=64QAM




%% Initialisation

SNR = zeros(1,length(SNRindB));

r_zf          = zeros(L,1);
r_zf_osic   = zeros(L,1);
r_mmse        = zeros(L,1);
r_mmse_osic = zeros(L,1);
   
Erreur_zf           = zeros(1,length(SNRindB));
Erreur_zf_osic    = zeros(1,length(SNRindB));
Erreur_mmse         = zeros(1,length(SNRindB));
Erreur_mmse_osic  = zeros(1,length(SNRindB));
 
for index = 1:length(SNRindB)
    SNR(index) = 10^(SNRindB(index)/10);
    sigma = sqrt(1/SNR(index));
    
    for simcnt= 1:Reps
        MsgBits  = randi([0,1], [log2(M),L]);   % Génération de bits à convertir en symboles
        Msg_symbole = bi2de(MsgBits', 'left-msb');  % Conversion des données binaires en symboles
        Msg_Mod = qammod(Msg_symbole, M);  % Modulation
        
        % Arrangement des symboles pour les transmettre un par un sur Nt antennes
        x = reshape(Msg_Mod,Nt,L/Nt);      
        
        %% Traitement à la réception
        H = complex(rand(Nr,Nt),rand(Nr,Nt));
        r = awgn(H*x,SNRindB(index),'measured')  ;  % Signal recu                                       
        for col_idx = 1:L/Nt

            rsic = r(:,col_idx);
            r_zf((col_idx-1)*Nt+1:col_idx*Nt)   = zf(rsic,H,M);
            r_zf_osic((col_idx-1)*Nt+1:col_idx*Nt) = osic_zf(rsic,H,M);
            r_mmse((col_idx-1)*Nt+1:col_idx*Nt) = mmse(rsic,H,sigma,M);
            r_mmse_osic((col_idx-1)*Nt+1:col_idx*Nt) = osic_mmse(rsic,H,sigma,M);
            
        end
        r_zf_bin = (de2bi(r_zf,'left-msb'));
        r_zf_sic_bin = de2bi(r_zf_osic,'left-msb');
        r_mmse_bin = de2bi(r_mmse,'left-msb');
        r_mmse_sic_bin = de2bi(r_mmse_osic,'left-msb');
        
        Erreur_zf(index)      = Erreur_zf(index) + biterr(r_zf_bin,MsgBits');
        Erreur_zf_osic(index)   = Erreur_zf_osic(index) + biterr(r_zf_sic_bin,MsgBits');
        Erreur_mmse(index)    = Erreur_mmse(index) + biterr(r_mmse_bin,MsgBits');
        Erreur_mmse_osic(index) = Erreur_mmse_osic(index) + biterr(r_mmse_sic_bin,MsgBits');
        
    end 
end % fin du boucle SNR

TotalBits = ((L*log2(M))*Reps);
 BER_zf   = Erreur_zf./TotalBits;
 BER_zf_osic  = Erreur_zf_osic./TotalBits;
 BER_mmse = Erreur_mmse./TotalBits;
 BER_mmse_osic = Erreur_mmse_osic./TotalBits;

%% Visualisation graphique des résultats

semilogy(SNRindB,BER_zf,'-ro');hold on;
semilogy(SNRindB,BER_mmse,'-b*');hold on;
semilogy(SNRindB,BER_zf_osic,'-kv');hold on;
semilogy(SNRindB,BER_mmse_osic,'-mpentagram');hold on;
xlabel('SNR (dB)');ylabel('BER');
legend('ZF','MMSE','ZF-OSIC','MMSE-OSIC');
grid on;
title('Performance VBLAST 4x4 , 64-QAM')