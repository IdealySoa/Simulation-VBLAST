%% Capacité VBLAST
close all;
clear;
clc;
%% Paramètre de la simulation
Nt = 2:2:8;
Nr = Nt;
SNRindB = 0:2:40;

%% Initialisation 
C = zeros( length(Nt), length(SNRindB));
SNR = zeros(1,length(SNRindB));
%% 
for index = 1: length(Nt)
    H = complex(rand(Nr(index),Nt(index)),rand(Nr(index),Nt(index)));
    for i = 1:length(SNRindB) 
        SNR(i) = 10^(SNRindB(i)/10);
        C(index,i) = log2(det(eye(Nt(index))+(SNR(i)/Nt(index))*(H*H')));
    end
    
end
%% 
figure 
plot(SNRindB,C,'-*');
xlabel('SNR (dB)');ylabel('Capacité (bits/s/Hz)');
grid on;
legend('Nt=Nr=2','Nt=Nr=4','Nt=Nr=6','Nt=Nr=8');
title('Capacité VBLAST')
