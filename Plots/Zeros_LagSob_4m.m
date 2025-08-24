% LAGUERRE SOBOLEV Computes and compares computed roots to values available
% in literature.
% This code generates Table 2 and associated errors
clearvars
close all
%% Define the problem and construct Gauss quadrature rule
m=10; g=1; s=1; nd=[m m]; same=1;
alpha = -1/2;  a0=alpha+1;
ab=r_laguerre(m,alpha); 
zw=gauss(m,ab);
xw=[zw(:,1) zw(:,1) zw(:,2) g*zw(:,2)];


n = 10;
%% Discretized Stieltjes procedure [Gautschi, Zhang 1995]
B=stieltjes_sob(m,s,nd,xw,a0,same);

% Modified Chebyshev - hier fehlen noch die Richtigen Momente, siehe S.32
% Gautschi 
gamma = g;
mom=zeros(2,2*n);
mom(1,1) = 2; mom(2,1)=2*gamma;
abm = r_laguerre(2*n-1);
B_cheb = chebyshev_sob(n,mom,abm);

%% New proposed methods
% These are based on matrix manipulation, so first the starting vector and
% Jordan matrix are formed. These generate the Krylov subspace

% Starting vector
w = zeros(2*m,1); w(2:2:2*m) = sqrt(zw(:,2));
% Jordan matrix
Z = zeros(2*m);
Z(1:2:end,1:2:end) = diag(zw(:,1));
Z(2:2:end,2:2:end) = diag(zw(:,1));
for k=2*m:-1:1%1:2*m
    if mod(k,2)==1
       Z(k,k+1) = sqrt(g);
    end
end


% Arnoldi iteration
[V,H] = Arnoldi(Z,w,m+1);
% Updating procedure
Hup = updating(Z,w,'PR');

%% plot computed zeros

% Annahme: A, B, C und D sind Ihre Matrizen (Sie müssen sie entsprechend definieren)
A = H; % Beispiel für A
B = Hup; % Beispiel für B
C = B; % Beispiel für C
D = B_cheb; % Beispiel für D

% Initialisierung von Zellen für die Matrizen A_i, B_i, C_i und D_i
A_i_cells = cell(1, 10);
B_i_cells = cell(1, 10);
C_i_cells = cell(1, 10);
D_i_cells = cell(1, 10);

% Schleife über i von 1 bis 10
for i = 1:10
    % Extrahiere die Submatrix A_i der Größe i x i
    A_i = A(1:i, 1:i);
    
    % Extrahiere die Submatrix B_i der Größe i x i
    B_i = B(1:i, 1:i);
    
    % Extrahiere die Submatrix C_i der Größe i x i
    C_i = C(1:i, 1:i);
    
    % Extrahiere die Submatrix D_i der Größe i x i
    D_i = D(1:i, 1:i);
    
    % Speichere die Submatrizen A_i, B_i, C_i und D_i in den Zellen
    A_i_cells{i} = A_i;
    B_i_cells{i} = B_i;
    C_i_cells{i} = C_i;
    D_i_cells{i} = D_i;
end

% Plotten der Eigenwerte der Matrizen A_i, B_i, C_i und D_i
figure;

for i = 1:10
    % Berechnung der Eigenwerte für A_i
    eigenvalues_A_i = eig(A_i_cells{i});
    
    % Berechnung der Eigenwerte für B_i
    eigenvalues_B_i = eig(B_i_cells{i});
    
    % Berechnung der Eigenwerte für C_i
    eigenvalues_C_i = eig(C_i_cells{i});                      
    
    % Berechnung der Eigenwerte für D_i
    eigenvalues_D_i = eig(D_i_cells{i});
    
    % Diskrete Schritte für die y-Achse
    y_values_A_i = i * ones(size(eigenvalues_A_i));
    y_values_B_i = i * ones(size(eigenvalues_B_i));
    y_values_C_i = i * ones(size(eigenvalues_C_i));
    y_values_D_i = i * ones(size(eigenvalues_D_i));
    
    % Plot der Eigenwerte von A_i in Rot
    scatter(real(eigenvalues_A_i), y_values_A_i, 'rx');
    hold on;
    
    % Plot der Eigenwerte von B_i in Blau
    scatter(real(eigenvalues_B_i), y_values_B_i, 'bo');
    
    % Plot der Eigenwerte von C_i in Grün
    scatter(real(eigenvalues_C_i), y_values_C_i, 'go');
    
    % Plot der Eigenwerte von D_i in Pink
    %scatter(real(eigenvalues_D_i), y_values_D_i, 'm+');

     % Berechnung und Ausgabe der Konditionszahl für jede Matrix
    cond_A_i = cond(A_i);
    cond_B_i = cond(B_i);
    cond_C_i = cond(C_i);
    cond_D_i = cond(D_i);
    
    fprintf('Konditionszahl von A_%d: %.4f\n', i, cond_A_i);
    fprintf('Konditionszahl von B_%d: %.4f\n', i, cond_B_i);
    fprintf('Konditionszahl von C_%d: %.4f\n', i, cond_C_i);
    fprintf('Konditionszahl von D_%d: %.4f\n', i, cond_D_i);
end

% Achsentitel und Beschriftungen
title(['Arnoldi Rot, AufteilungS Blau, Stieltjes in Grün, Chebyshev in Pink']);
xlabel('Reale Achse');
ylabel('i');

% Gitter hinzufügen
grid off;

% Achsen anpassen
axis tight;

hold off;