clc ; 
clear ; 
close ;

f0 = 1e3; % Frequence du signal en Hz
fe = 1e5; % Frequence d'echantillonage (en Hz) >= 2 * fmax (f0)
T = 100/f0;
t = 0:1/fe:T; % Axe de temps : la durée du signal 
SNR = 10; % en dB
snr = 10^(SNR/10); % Valeur Linéaire 
N = 2; % Nombre du coefficient du filtre de Wiener h = [h(0), .... h(N-1)]

sk = sin(2*pi*f0*t(1:end-1)); % Signal Utile | t(1:end-1) => Nombre Paire
P_sk = var(sk); % puissance du signal sk
% La commande var calcul la variance => la puissance du signal 
sigma_n = sqrt(P_sk/snr); % Ecart type du bruit 
nk = sigma_n * randn(1,length(sk)); %Bruit de variance sigma_n 
rk = sk + nk; % Observation bruitée r(k)

%===================> Estimation de R_r et r_sr : 
R_r = zeros(N,N);
R_s = zeros(N,N);
R_n = zeros(N,N);
r_sr = zeros(N,1);

for iteration = 1:length(rk)-N+1
    R_r = R_r + rk(iteration:iteration+N-1)' * rk(iteration:iteration+N-1);
    R_s = R_s + sk(iteration:iteration+N-1)' * sk(iteration:iteration+N-1);
    R_n = R_n + nk(iteration:iteration+N-1)' * nk(iteration:iteration+N-1);
    r_sr = r_sr + sk(iteration) * rk(iteration:iteration+N-1)'; 
end 
R_r_est = R_r / length(rk)-N+1; % length(rk)-N+1 = M (Ecrites dans les notes)
r_sr_est = r_sr / length(rk)-N+1;
R_s_est = R_s / length(rk)-N+1;
R_n_est = R_n / length(rk)-N+1;
h_opt = inv(R_r_est)*r_sr_est; 
d_hat = filter(h_opt,1,rk); % RIF : produit Convuktionnel entre rk et rif


%====================> Affichage : 
figure(1);     
plot(t(1:end-1),sk,'-k','LineWidth',2)
hold on 
plot(t(1:end-1),rk,'-r','LineWidth',2)
hold on 
plot(t(1:end-1),d_hat,'-b','LineWidth',2)
xlabel('Temps [Sec]'), ylabel('Amplitude')
legend('s(k)','r(k)')
axis([0 t(end) -2 2]), grid 
% Zoom pour voir plus de details | On essel avec SNR << : bruit moins
% variable 


SNR_in = SNR; % On fixe l'entrée 
SNR_out = 10*log10(var(filter(h_opt,1,sk)))/var(filter(h_opt,1,nk)); % On calcul la sortie 