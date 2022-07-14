%Program: Détecteur de Bayes (Perr, Pfa, Pnd) Vs SNR
%Date : 29/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% TP 1 : Detecteur de BAYES %%%%%%%%%%%%%%%%%%%%%%%%
clear all ; 
close all ; 
clc ; 

%========= Paramétres de simulation : 
N_tirage = 1e4; % Nombre de symbole s0 ou s1 par segment 
SNR = [0:13]; %Signal to noire Ratio en dB
snr = 10.^(SNR/10); % SNR en valeur linéaire
s0 = -1; %Signal pour bit transmis 0 
s1 =  1; %Signal pour bit transmis 1
p0 = 1/5 ; % Probabilité de signal émis s0 
p1 = 4/5 ; % Probabilité de signal émis s1
p_soi = 1; % Puissance du signal utile
sigma_n = sqrt(p_soi./snr); % Ecart Type du bruit 

%========= Detecteur de Bayes Pour un seul Echantillon : 
% Initialisation : 
Perr_hat = zeros(1,length(SNR));
Perr_hat_bis = zeros(1,length(SNR));
Pfa_hat = zeros(1,length(SNR));
Pnd_hat = zeros(1,length(SNR));
x_seuil = zeros(1,length(SNR));
Perr = zeros(1,length(SNR));

for i_SNR =1:length(SNR)
    Nb_segement = 0; % Nombre de segment Tx pour claculer Perr
    Nb_segement_FA = 0; %Npmbre de Segment Tx pour calculer Pfa
    Nb_segement_ND = 0; % Nombre de Segment Tx pour calculer Pnd
    
    Nb_erreur = 0; % Nombre d'erreur total pour calculer Perr
    Nb_erreur_ND = 0; % Nombre d'erreur total pour calculer Pnd
    Nb_erreur_FA = 0; % Nombre d'erreur total pour calculer Pfa
    
    temp_NB_erreur = 0; %Nomre d'erreurs par segment 
    
    %%%%%%%% Perr %%%%%%%%
    while Nb_erreur < 300  
        Vecteur_s = (-1).^random('Binomial',1,p0,[1,N_tirage]);% Génération du signal Utile 
        Vecteur_n = sigma_n(i_SNR)*randn(1,N_tirage); %Génération de bruit de var sigma_n^2
        Vecteur_x = Vecteur_s + Vecteur_n;
        
        x_seuil(i_SNR) = (((s0+s1)/2)+(log(p0/p1) / (s1-s0)) * sigma_n(i_SNR)^2); % Seuil de décision 
        Vecteur_s_hat = (-1).^double(Vecteur_x < x_seuil(i_SNR)); % Detection du symbole transmis 
        
        temp_Nb_erreur = length(Vecteur_s) - sum(Vecteur_s_hat == Vecteur_s); % Calcul du nombre d'erreurs par segement
                     
        Nb_erreur = Nb_erreur + temp_Nb_erreur; % nombre d'erreur cumul�
        Nb_segement = Nb_segement + 1; %Compteur de segements  
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%> Pfa <%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while Nb_erreur_FA < 500
        
        Vecteur_s = (-1).^random('Binomial',1,p0,[1, N_tirage]); % G�n�ration du signal utile
        Vecteur_n = sigma_n(i_SNR) * randn(1,N_tirage); %% G�n�ration du bruit de variance sigma_n^2
        Vecteur_x = Vecteur_s + Vecteur_n;% G�n�ration du signal utile bruit� de longueur N_tirage = 10^4

        x_seuil(i_SNR) = (((s0+s1)/2) + (log(p0/p1) / (s1-s0)) * sigma_n(i_SNR)^2);% Seuil de d�cision
        Vecteur_s_hat =  (-1).^double(Vecteur_x < x_seuil(i_SNR));% D�tection du symbole transmis s0 ou s1
        
        Position_Erreur_FA = double(Vecteur_s_hat == Vecteur_s);% Vecteur indicateur d'erreurs f.a quand = 0

        
        temp_Nb_erreur_FA = 0;
        for i_erreur = 1:length(Position_Erreur_FA)
           if (Position_Erreur_FA(i_erreur) == 0) && (Vecteur_s(i_erreur) == -1)% Condition pour une erreur de type -1 --> +1
               temp_Nb_erreur_FA = temp_Nb_erreur_FA + 1;%Compteur d'erreurs de type -1 --> +1 par segements
           end
        end
        
        
        Nb_erreur_FA = Nb_erreur_FA + temp_Nb_erreur_FA; %Compteur d'erreurs de type -1 --> +1 total
        Nb_segement_FA = Nb_segement_FA + 1;%Compteur de segements
        
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%> Pnd <%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while Nb_erreur_ND < 500
        
        Vecteur_s = (-1).^random('Binomial',1,p0,[1, N_tirage]); % G�n�ration du signal utile
        Vecteur_n = sigma_n(i_SNR) * randn(1,N_tirage); % G�n�ration du bruit de variance sigma_n^2
        Vecteur_x = Vecteur_s + Vecteur_n;%% G�n�ration du signal utile bruit� de longueur N_tirage = 10^4

        x_seuil(i_SNR) = (((s0+s1)/2) + (log(p0/p1) / (s1-s0)) * sigma_n(i_SNR)^2);% Seuil de d�cision
        Vecteur_s_hat =  (-1).^double(Vecteur_x < x_seuil(i_SNR));% D�tection du symbole transmis s0 ou s1
        
        Position_Erreur_ND = double(Vecteur_s_hat == Vecteur_s);% Vecteur indicateur d'erreurs n.d quand = 0
        
        temp_Nb_erreur_ND = 0;
        for i_erreur = 1:length(Position_Erreur_ND)
           if (Position_Erreur_ND(i_erreur) == 0) && (Vecteur_s(i_erreur) == 1)% Condition pour une erreur de type +1 --> -1
               temp_Nb_erreur_ND = temp_Nb_erreur_ND + 1;%Compteur d'erreurs de type +1 --> -1 par segements
           end
        end
        
        
        Nb_erreur_ND = Nb_erreur_ND + temp_Nb_erreur_ND;%Compteur d'erreurs de type +1 --> -1 total
        Nb_segement_ND = Nb_segement_ND + 1;%Compteur de segements
        
    end
    
    Perr_hat(i_SNR) = (Nb_erreur) / (Nb_segement * N_tirage);
    
    Pfa_hat(i_SNR) = (Nb_erreur_FA) / (Nb_segement_FA * N_tirage) / p0; % Attention ! il faut diviser par p0... voir la condition "if" plus haut
    Pnd_hat(i_SNR) = (Nb_erreur_ND) / (Nb_segement_ND * N_tirage) / p1; % Attention ! il faut diviser par p1... voir la condition "if" plus haut
    
    
    Perr(i_SNR) = p0 * qfunc((x_seuil(i_SNR)-s0)/sigma_n(i_SNR)) + p1 * qfunc((s1-x_seuil(i_SNR))/sigma_n(i_SNR)); % Perr th�orique
    Pnd(i_SNR) = qfunc((s1-x_seuil(i_SNR))/sigma_n(i_SNR)); % Pnd th�orique
    Pfa(i_SNR) = qfunc((x_seuil(i_SNR)-s0)/sigma_n(i_SNR)); % Pfa th�orique
    
%     waitbar(i_SNR/length(SNR));
    
end


figure(1)
semilogy(SNR,Perr_hat,'*','LineWidth',2)
hold on
semilogy(SNR,Perr,'LineWidth',2)


xlabel('SNR [dB]'),ylabel('Probabilit� d''erreur')
legend('P_{err} empirique','P_{err} th�orique')
axis([0 13 1e-6 1e0]),grid


ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(2)
semilogy(SNR,Pnd_hat,'*','LineWidth',2)
hold on
semilogy(SNR,Pnd,'LineWidth',2)


xlabel('SNR [dB]'),ylabel('Probabilit�')
legend('P_{nd} empirique','P_{nd} th�orique')
axis([0 13 1e-6 1e0]),grid


ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(3)
semilogy(SNR,1 - Pnd_hat,'*','LineWidth',2)
hold on
semilogy(SNR,1 - Pnd,'LineWidth',2)


xlabel('SNR [dB]'),ylabel('Probabilit�')
legend('P_{d} empirique','P_{d} th�orique')
axis([0 13 1e-6 1.2e0]),grid


ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';


figure(4)
semilogy(SNR,Pfa_hat,'*','LineWidth',2)
hold on
semilogy(SNR,Pfa,'LineWidth',2)


xlabel('SNR [dB]'),ylabel('Probabilit�')
legend('P_{fa} empirique','P_{fa} th�orique')
axis([0 13 1e-6 1e0]),grid


ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';



    
