%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Program: Détecteur de Bayes avec un seul échantillon par symbole
%Date : 29/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
%===========> Paramètres de simulation
N_tirage = 1e4;% Nombre de symbole s0 ou s1 transmis par segment
SNR = [0,1,2,3,4,5]; % [dB]
snr = 10.^(SNR/10); % SNR en Valeur linéaire
s0 = -1;% signal pour le bit '0'
s1 = +1;% signal pour le bit '1'
p0 = 1/2; % Probabilité s_0
p1 = 1/2; % Probabilité s_1
p_soi = 1;% puissance signal utile
sigma_n = sqrt(p_soi./snr); % écart type du bruit
Vec_x_seuil = [-1:0.1:1];% seuil de décision variant de -1 à 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Détecteur de Bayes à un seul échantillon pour chaque décision       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================> Initialisation <===============
Pfa_hat  = zeros(length(Vec_x_seuil),length(SNR));
Pd_hat = zeros(length(Vec_x_seuil),length(SNR));

Pfa  = zeros(length(Vec_x_seuil),length(SNR));
Pd = zeros(length(Vec_x_seuil),length(SNR));
%==================================================




for i_SNR = 1:length(SNR)
    
    for i_seuil = 1:length(Vec_x_seuil)
        
       
        Nb_segement_FA = 0;% Nombre de segments TX pour calculer Pfa
        Nb_segement_ND = 0;% Nombre de segments TX pour calculer Pnd

        
        Nb_erreur_ND = 0; % Nombre d'erreurs total pour calculer Pnd
        Nb_erreur_FA = 0; % Nombre d'erreurs total pour calculer Pfa
        temp_Nb_erreur = 0; % Nombre d'erreurs par segment

        %%%%%%%%%%%%%%%%%%%%%%%%%%%> Pfa <%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        while Nb_erreur_FA < 500

            Vecteur_s = (-1).^random('Binomial',1,p0,[1, N_tirage]); % Génération du signal utile
            Vecteur_n = sigma_n(i_SNR) * randn(1,N_tirage); % Génération du bruit de variance sigma_n^2
            Vecteur_x = Vecteur_s + Vecteur_n;% % Génération du signal utile bruité de longueur N_tirage = 10^4

            Vecteur_s_hat =  (-1).^double(Vecteur_x < Vec_x_seuil(i_seuil)); % Détection du symbole transmis s0 ou s1

            Position_Erreur_FA = double(Vecteur_s_hat == Vecteur_s); % Vecteur indicateur d'erreurs f.a quand = 0

            temp_Nb_erreur_FA = 0;
            for i_erreur = 1:length(Position_Erreur_FA)
               if (Position_Erreur_FA(i_erreur) == 0) && (Vecteur_s(i_erreur) == -1) % Condition pour une erreur de type -1 --> +1
                   temp_Nb_erreur_FA = temp_Nb_erreur_FA + 1;
               end
            end

            Nb_erreur_FA = Nb_erreur_FA + temp_Nb_erreur_FA;
            Nb_segement_FA = Nb_segement_FA + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%> Pd <%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        while Nb_erreur_ND < 500

            Vecteur_s = (-1).^random('Binomial',1,p0,[1, N_tirage]); % Génération du signal utile
            Vecteur_n = sigma_n(i_SNR) * randn(1,N_tirage); % Génération du bruit de variance sigma_n^2
            Vecteur_x = Vecteur_s + Vecteur_n;% Génération du signal utile bruité de longueur N_tirage = 10^4

            Vecteur_s_hat =  (-1).^double(Vecteur_x < Vec_x_seuil(i_seuil)); % Détection du symbole transmis s0 ou s1

            Position_Erreur_ND = double(Vecteur_s_hat == Vecteur_s);% Vecteur indicateur d'erreurs n.d quand = 0

            temp_Nb_erreur_ND = 0;
            for i_erreur = 1:length(Position_Erreur_ND)
               if (Position_Erreur_ND(i_erreur) == 0) && (Vecteur_s(i_erreur) == 1) % Condition pour une erreur de type +1 --> -1
                   temp_Nb_erreur_ND = temp_Nb_erreur_ND + 1;
               end
            end


            Nb_erreur_ND = Nb_erreur_ND + temp_Nb_erreur_ND;
            Nb_segement_ND = Nb_segement_ND + 1;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    Pfa_hat(i_seuil,i_SNR) = (Nb_erreur_FA) / (Nb_segement_FA * N_tirage) / p0; % Attention ! il faut diviser par p0... voir la condition "if" plus haut
    Pd_hat(i_seuil,i_SNR) = 1 -  ((Nb_erreur_ND) / (Nb_segement_ND * N_tirage) / p1); % Attention ! il faut diviser par p1... voir la condition "if" plus haut
    
    Pd(i_seuil, i_SNR) = 1 - qfunc((s1-Vec_x_seuil(i_seuil))/sigma_n(i_SNR)); % Pd théorique
    Pfa(i_seuil, i_SNR) = qfunc((Vec_x_seuil(i_seuil)-s0)/sigma_n(i_SNR));% % Pfa théorique
    
    end
    
    
    
end


figure(1)
plot(Vec_x_seuil,Pd_hat,'LineWidth',2)
xlabel('x_{seuil}'),ylabel('P_{d}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([-1 1 0.5 1e0]),grid
title('Courbes empiriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(11)
plot(Vec_x_seuil,Pd,'LineWidth',2)
xlabel('x_{seuil}'),ylabel('P_{d}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([-1 1 0.5 1e0]),grid
title('Courbes théoriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';


figure(2)
plot(Vec_x_seuil,Pfa_hat,'LineWidth',2)
xlabel('x_{seuil}'),ylabel('P_{fa}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([-1 1 0 0.5]),grid
title('Courbes empiriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(22)
plot(Vec_x_seuil,Pfa,'LineWidth',2)
xlabel('x_{seuil}'),ylabel('P_{fa}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([-1 1 0 0.5]),grid
title('Courbes théoriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(3)
plot(Pfa_hat,Pd_hat,'LineWidth',2)
xlabel('P_{fa}'),ylabel('P_{d}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([0 0.5 0.5 1]),grid
title('Courbes COR empiriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';

figure(33)
plot(Pfa,Pd,'LineWidth',2)
xlabel('P_{fa}'),ylabel('P_{d}')
legend('SNR = 0dB','SNR = 1dB','SNR = 2dB','SNR = 3dB','SNR = 4dB','SNR = 5dB')
axis([0 0.5 0.5 1e0]),grid
title('Courbes COR théoriques')
ax = gca;
ax.FontSize = 16;
ax.TitleFontWeight = 'bold';
