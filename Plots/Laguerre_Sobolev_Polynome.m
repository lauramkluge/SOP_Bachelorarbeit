% Laguerre-Sobolev Polynome 
close all
clearvars
n = 10; s = 1; nd = [n n]; same = 1;
alpha = 2;
a0 = alpha + 1;
ab = r_laguerre(n,alpha); 
iter = 10;

%gamma_values = 10.^(-6:10);
gamma_values = [10^4];

% Initialisiere Array für maximale EW-Konditionszahlen von ...
max_cond_Arn_iter = zeros(1, iter); % Arnoldi
max_cond_Auf_iter = zeros(1, iter); % Aufteilungsverfahren
max_cond_Sti_iter = zeros(1, iter); % Stieltjes-Verfahren



for gamma_index = 1:length(gamma_values)
    gamma = gamma_values(gamma_index);
    zw = gauss(n, ab);
    
    %% Stieltjes-Verfahren von Gautschi und Zhang
    xw = [zw(:, 1) zw(:, 1) zw(:, 2) gamma * zw(:, 2)];
    [Stieltjes, Norm_Stieltjes] = stieltjes_sob(n, s, nd, xw, a0, same);

    %% Verfahren von Buggenhout
    % Gewichtsvektor und Jordan-ähnliche Matrix berechnen
    w = zeros(2 * n, 1);
    w(2:2:2 * n) = sqrt(zw(:, 2));
    A = zeros(2 * n);
    A(1:2:end, 1:2:end) = diag(zw(:, 1));
    A(2:2:end, 2:2:end) = diag(zw(:, 1));
    for k = 2 * n:-1:1
        if mod(k, 2) == 1
            A(k, k + 1) = sqrt(gamma);
        end
    end

    % Arnoldi-Verfahren
    [V, H_Arnoldi] = Arnoldi(A, w, n + 1);
    % Aufteilungsverfahren
    H_Aufteilung = updating(A, w, 'PR');

    %% Abbildungen erstellen

    % Initialisiere Abweichungsmatrizen
    ArnSti_Abw = zeros(1,iter);
    AufSti_Abw = zeros(1,iter);
    ArnAuf_Abw = zeros(1,iter);
    
    figure('Name',"Nullstellen der orthonormalen Polynome p_k für gamma="+num2str(gamma),'NumberTitle','off');

    for k = 1:iter
        % Berechnung der Eigenwerte und Speicherung der Konditionszahl X_k
        [X_Arnoldi, D_Arnoldi] = eig(H_Arnoldi(1:k, 1:k), "vector");
        [X_Aufteilung, D_Aufteilung] = eig(H_Aufteilung(1:k, 1:k), "vector");
        [X_Stieltjes, D_Stieltjes] = sobzeros_mine(k,k,Stieltjes, sqrt(Norm_Stieltjes(1:k)));

        max_cond_Arn_iter(k) = cond(X_Arnoldi);
        max_cond_Auf_iter(k) = cond(X_Aufteilung);
        max_cond_Sti_iter(k) = cond(X_Stieltjes);

        D_Arnoldi = sort(D_Arnoldi, 'ComparisonMethod','real');
        D_Aufteilung = sort(D_Aufteilung, 'ComparisonMethod','real'); % bei HH können imaginäre Anteile entstehen 
        D_Stieltjes = sort(D_Stieltjes, 'ComparisonMethod','real');

        % Plot der Lage der Nullstellen
        scatter(real(D_Arnoldi), k * ones(size(D_Arnoldi)), 'rx');
        hold on;
        scatter(real(D_Aufteilung), k * ones(size(D_Aufteilung)), 'b+');
        scatter(real(D_Stieltjes), k * ones(size(D_Stieltjes)), 'go');
        
        % Abweichungen der Ergebnisse zwischen den Verfahren
        ArnSti_Abw(k) = max(abs(D_Arnoldi - D_Stieltjes));
        AufSti_Abw(k) = max(abs(D_Aufteilung - D_Stieltjes));
        ArnAuf_Abw(k) = max(abs(D_Arnoldi - D_Aufteilung));
    end

    % % Nullstellen-Plot
    set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14)
    xlabel('Reelle Achse');
    ylabel('Polynomgrad');
    % lgd = legend('Arnoldi-', 'Aufteilungs-', 'Stieltjes-');
    % title(lgd,'Verfahren');
    grid off;
    axis tight;
    hold off;

    % Abweichungs-Plot
    figure('Name',"Nullstellenabweichung für p_k für gamma = "+num2str(gamma),'NumberTitle','off');    
    semilogy(AufSti_Abw, '-', 'LineWidth', 3);
    hold on;
    semilogy(ArnSti_Abw, ':', 'LineWidth', 3);
    semilogy(ArnAuf_Abw, '-.', 'LineWidth', 3);
    set(gca, 'FontName', 'Times New Roman', 'Fontsize', 18)
    xlabel('k');
    ylabel('Abweichung')
    % legend('Arnoldi/Stieltjes', 'Aufteilung/Stieltjes', 'Arnoldi/Aufteilung');
    hold off;

    % Konditionszahl-Plot
    figure('Name',"Konditionszahlen der Eigenwertberechnungen von p_k für gamma = "+num2str(gamma),'NumberTitle','off');
    semilogy(max_cond_Arn_iter, 'r:', 'LineWidth', 3);
    hold on;
    semilogy(max_cond_Auf_iter, 'b--', 'LineWidth', 3);
    semilogy(max_cond_Sti_iter, 'g-.', 'LineWidth', 3);
    set(gca, 'FontName', 'Times New Roman', 'Fontsize', 18)
    xlabel('k');
    ylabel('Konditionszahl');
    % legend('Arnoldi-', 'Aufteilungs-', 'Stieltjes');
    % title(lgd,'Verfahren');
    hold off;
end