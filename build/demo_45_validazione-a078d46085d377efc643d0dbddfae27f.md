# Demo 45 - Validazione di un modello VAR

## Scopo

La demo simula un VAR(5), seleziona l'ordine con AIC e confronta un modello sottostimato con quello selezionato. Il messaggio centrale è che la stabilità non basta: un modello può essere stabile e lasciare struttura temporale nei residui.

## Risultati attesi

- l'AIC seleziona un ordine vicino a quello reale, con possibili variazioni dovute al campionamento;
- il VAR(2) mostra autocorrelazione residua maggiore e una consistenza peggiore;
- entrambi i modelli possono risultare stabili, ma soltanto il controllo congiunto dei residui consente di giudicare l'adeguatezza.

I valori numerici cambiano con le realizzazioni casuali. Lo script imposta il generatore quando è necessario per rendere il confronto riproducibile.

## Output verificato

![Output MATLAB di Demo 45 - Validazione di un modello VAR](../assets/demo/demo_45_01.png)

![Output MATLAB di Demo 45 - Validazione di un modello VAR](../assets/demo/demo_45_02.png)

## Esecuzione

Scaricare [Demo_45_validazione.m](../matlab/Demo_45_validazione.m), aprirlo in MATLAB ed eseguirlo per sezioni oppure integralmente. Gli script sono stati verificati con MATLAB R2025b.

## Codice MATLAB

```matlab
%% DEMO 45 - VALIDAZIONE DI UN MODELLO VAR MULTIVARIATO
% La demo genera un sistema VAR(5) a quattro variabili, stima l'ordine con
% il criterio AIC e confronta due modelli:
%   1) un VAR(2), con ordine volutamente sottostimato;
%   2) un VAR(p_AIC), con ordine selezionato dai dati.
%
% Per entrambi i modelli vengono valutati:
%   - bianchezza dei residui: test portmanteau multivariato di Ljung-Box;
%   - stabilita': autovalori della matrice companion;
%   - consistenza: confronto tra le correlazioni dei dati reali e simulati.
%
% Convenzione: A(i,j,k) descrive l'influenza j -> i al ritardo k.
% Tutto il codice e' contenuto in questo script, senza funzioni locali.

clearvars;
close all;
clc;

rng(12, "twister");

%% 1. PARAMETRI DELLA SIMULAZIONE
M = 4;                         % Numero di variabili
p_true = 5;                    % Ordine reale del processo
N = 2000;                      % Campioni mantenuti
burn_in = 500;                 % Campioni eliminati come transitorio
alpha = 0.05;                  % Livello di significativita'
h_test = 20;                   % Lag massimo per portmanteau e consistenza

% Covarianza delle innovazioni.
Sigma_u = [1.0  0.2  0.0  0.0;
           0.2  0.8  0.1  0.0;
           0.0  0.1  1.2  0.15;
           0.0  0.0  0.15 0.9];

%% 2. DEFINIZIONE DEL VAR(5) REALE
% x(t) = A1*x(t-1) + ... + A5*x(t-5) + u(t)
A_true = zeros(M, M, p_true);

% Dinamiche al primo ritardo.
A_true(:,:,1) = [0.55 0.00 0.00 0.00;
                 0.35 0.45 0.00 0.00;
                 0.00 0.25 0.40 0.00;
                 0.00 0.00 0.00 0.50];

% Componenti autoregressive al secondo ritardo.
A_true(:,:,2) = diag([-0.18 -0.12 -0.15 -0.16]);

% Connessioni presenti soltanto al quinto ritardo.
% Questi termini non possono essere descritti correttamente da un VAR(2).
A_true(4,1,5) =  0.45;         % x1 -> x4
A_true(3,4,5) = -0.30;         % x4 -> x3
A_true(2,2,5) =  0.20;         % memoria di x2 al quinto ritardo

% Verifica della stabilita' del modello generatore.
A_true_top = reshape(A_true, M, M*p_true);
A_true_bottom = [eye(M*(p_true-1)), zeros(M*(p_true-1), M)];
A_true_companion = [A_true_top; A_true_bottom];
rho_true = max(abs(eig(A_true_companion)));

fprintf('Raggio spettrale del VAR(5) reale: %.4f\n', rho_true);
if rho_true >= 1
    error('Il modello generatore non e'' stabile.');
end

%% 3. GENERAZIONE DEI DATI SINTETICI
n_total = N + burn_in + p_true;
U = chol(Sigma_u, "lower") * randn(M, n_total);
X_all = zeros(M, n_total);
X_all(:,1:p_true) = randn(M, p_true);

for t = p_true+1:n_total
    X_all(:,t) = U(:,t);
    for lag = 1:p_true
        X_all(:,t) = X_all(:,t) ...
            + A_true(:,:,lag) * X_all(:,t-lag);
    end
end

X = X_all(:, burn_in+p_true+1:end);

%% 4. SELEZIONE DELL'ORDINE MEDIANTE AIC
% Per ogni ordine candidato viene stimato un VAR e calcolato:
%
% AIC(p) = log(det(Sigma_u_hat)) + 2*k(p)/(N-p)
%
% dove k(p) = M^2*p + M e' il numero totale di coefficienti stimati,
% includendo le M intercette. Si seleziona l'ordine che minimizza l'AIC.

candidate_orders = 1:10;
AIC_values = zeros(size(candidate_orders));

for candidate_index = 1:length(candidate_orders)
    p_candidate = candidate_orders(candidate_index);
    n_candidate = N - p_candidate;

    Y_candidate = X(:, p_candidate+1:N);
    Z_candidate = zeros(1 + M*p_candidate, n_candidate);
    Z_candidate(1,:) = 1;

    for lag = 1:p_candidate
        rows = 2 + (lag-1)*M : 1 + lag*M;
        Z_candidate(rows,:) = X(:, p_candidate+1-lag:N-lag);
    end

    B_candidate = Y_candidate / Z_candidate;
    residuals_candidate = Y_candidate - B_candidate * Z_candidate;
    Sigma_candidate = (residuals_candidate * residuals_candidate') ...
                    / n_candidate;

    % log(det(.)) viene calcolato dagli autovalori per evitare instabilita'
    % numeriche associate al determinante diretto.
    Sigma_candidate = (Sigma_candidate + Sigma_candidate') / 2;
    eigenvalues_candidate = max(eig(Sigma_candidate), realmin);
    log_determinant = sum(log(eigenvalues_candidate));

    number_parameters = M^2 * p_candidate + M;
    AIC_values(candidate_index) = log_determinant ...
        + 2 * number_parameters / n_candidate;
end

[~, minimum_AIC_index] = min(AIC_values);
p_AIC = candidate_orders(minimum_AIC_index);

fprintf('\nOrdine reale del processo:       %d\n', p_true);
fprintf('Ordine selezionato mediante AIC: %d\n', p_AIC);

figure("Color", "w", "Name", "Selezione dell'ordine mediante AIC");
plot(candidate_orders, AIC_values, "o-", "LineWidth", 1.4, ...
    "MarkerFaceColor", [0.2 0.45 0.75]);
hold on;
plot(p_AIC, AIC_values(minimum_AIC_index), "rp", ...
    "MarkerSize", 13, "MarkerFaceColor", "r");
grid on;
xticks(candidate_orders);
xlabel("Ordine candidato p");
ylabel("AIC(p)");
title(sprintf("Ordine selezionato dall'AIC: p = %d", p_AIC));

%% 5. STIMA E VALIDAZIONE CON DUE ORDINI DIVERSI
% Il primo ordine e' volutamente troppo basso; il secondo e' selezionato
% automaticamente mediante AIC.
model_orders = [2 p_AIC];
n_models = length(model_orders);

Q_portmanteau = zeros(1, n_models);
df_portmanteau = zeros(1, n_models);
p_portmanteau = zeros(1, n_models);
spectral_radius = zeros(1, n_models);
stability_index = zeros(1, n_models);
percent_consistency = zeros(1, n_models);
residual_correlation_strength = zeros(n_models, h_test);
residuals_models = cell(1, n_models);

for model_index = 1:n_models

    p_est = model_orders(model_index);
    n_effective = N - p_est;

    % ---------------------------------------------------------------------
    % 5.1 Stima multivariata ai minimi quadrati
    % ---------------------------------------------------------------------
    Y = X(:, p_est+1:N);
    Z = zeros(1 + M*p_est, n_effective);
    Z(1,:) = 1;

    for lag = 1:p_est
        rows = 2 + (lag-1)*M : 1 + lag*M;
        Z(rows,:) = X(:, p_est+1-lag:N-lag);
    end

    % Y = B*Z + E, con B = [intercetta A1 ... Ap].
    B_hat = Y / Z;
    residuals = Y - B_hat * Z;
    residuals_models{model_index} = residuals;

    A_hat = zeros(M, M, p_est);
    for lag = 1:p_est
        columns = 2 + (lag-1)*M : 1 + lag*M;
        A_hat(:,:,lag) = B_hat(:,columns);
    end

    % ---------------------------------------------------------------------
    % 5.2 Test portmanteau multivariato di Ljung-Box
    % ---------------------------------------------------------------------
    % H0: C_1 = ... = C_h = 0.
    % Il test considera congiuntamente auto- e cross-correlazioni residue.
    n_res = size(residuals, 2);
    C0 = (residuals * residuals') / n_res;
    C0_inv = pinv(C0);
    residual_std = sqrt(max(diag(C0), eps));
    residual_normalization = residual_std * residual_std';

    Q_value = 0;
    for lag = 1:h_test
        C_lag = (residuals(:,lag+1:end) ...
               * residuals(:,1:end-lag)') / n_res;

        Q_value = Q_value ...
            + trace(C_lag' * C0_inv * C_lag * C0_inv) / (n_res-lag);

        % Quantita' descrittiva usata nel grafico finale.
        R_lag = C_lag ./ residual_normalization;
        residual_correlation_strength(model_index,lag) = norm(R_lag, "fro");
    end

    Q_value = n_res * (n_res + 2) * Q_value;
    degrees_freedom = M^2 * (h_test - p_est);
    p_value = 1 - chi2cdf(Q_value, degrees_freedom);

    Q_portmanteau(model_index) = Q_value;
    df_portmanteau(model_index) = degrees_freedom;
    p_portmanteau(model_index) = p_value;

    % ---------------------------------------------------------------------
    % 5.3 Stabilita' del modello stimato
    % ---------------------------------------------------------------------
    if p_est == 1
        A_companion = A_hat(:,:,1);
    else
        A_top = reshape(A_hat, M, M*p_est);
        A_bottom = [eye(M*(p_est-1)), zeros(M*(p_est-1), M)];
        A_companion = [A_top; A_bottom];
    end

    spectral_radius(model_index) = max(abs(eig(A_companion)));
    stability_index(model_index) = log(spectral_radius(model_index));

    % ---------------------------------------------------------------------
    % 5.4 Percent consistency
    % ---------------------------------------------------------------------
    if spectral_radius(model_index) < 1
        Sigma_u_hat = (residuals * residuals') / n_res;
        Sigma_u_hat = (Sigma_u_hat + Sigma_u_hat') / 2;

        % Fattorizzazione robusta della covarianza stimata.
        [V_sigma, D_sigma] = eig(Sigma_u_hat);
        eigenvalues_sigma = max(diag(D_sigma), 0);
        innovation_factor = V_sigma * diag(sqrt(eigenvalues_sigma));

        n_sim_total = N + burn_in;
        U_sim = innovation_factor * randn(M, n_sim_total);
        X_sim_all = zeros(M, n_sim_total);
        X_sim_all(:,1:p_est) = X(:,1:p_est);
        intercept_hat = B_hat(:,1);

        for t_sim = p_est+1:n_sim_total
            X_sim_all(:,t_sim) = intercept_hat + U_sim(:,t_sim);
            for lag = 1:p_est
                X_sim_all(:,t_sim) = X_sim_all(:,t_sim) ...
                    + A_hat(:,:,lag) * X_sim_all(:,t_sim-lag);
            end
        end

        X_sim = X_sim_all(:, burn_in+1:burn_in+N);

        X_real_centered = X - mean(X, 2);
        X_sim_centered = X_sim - mean(X_sim, 2);

        C0_real = (X_real_centered * X_real_centered') / N;
        C0_sim = (X_sim_centered * X_sim_centered') / N;

        std_real = sqrt(max(diag(C0_real), eps));
        std_sim = sqrt(max(diag(C0_sim), eps));
        normalization_real = std_real * std_real';
        normalization_sim = std_sim * std_sim';

        R_real = zeros(M^2, h_test+1);
        R_sim = zeros(M^2, h_test+1);

        for lag = 0:h_test
            C_real_lag = (X_real_centered(:,lag+1:end) ...
                        * X_real_centered(:,1:end-lag)') / (N-lag);
            C_sim_lag = (X_sim_centered(:,lag+1:end) ...
                       * X_sim_centered(:,1:end-lag)') / (N-lag);

            R_real_lag = C_real_lag ./ normalization_real;
            R_sim_lag = C_sim_lag ./ normalization_sim;

            R_real(:,lag+1) = R_real_lag(:);
            R_sim(:,lag+1) = R_sim_lag(:);
        end

        percent_consistency(model_index) = 100 * ...
            (1 - norm(R_sim(:) - R_real(:), 2) / norm(R_real(:), 2));
    else
        percent_consistency(model_index) = NaN;
    end
end

%% 6. RISULTATI NUMERICI
portmanteau_result = strings(n_models,1);
stability_result = strings(n_models,1);

for model_index = 1:n_models
    if p_portmanteau(model_index) < alpha
        portmanteau_result(model_index) = "Residui non bianchi";
    else
        portmanteau_result(model_index) = "Residui compatibili con bianchezza";
    end

    if spectral_radius(model_index) < 1
        stability_result(model_index) = "Stabile";
    else
        stability_result(model_index) = "Non stabile";
    end
end

results = table(model_orders', Q_portmanteau', df_portmanteau', ...
    p_portmanteau', portmanteau_result, spectral_radius', ...
    stability_index', stability_result, percent_consistency', ...
    'VariableNames', {'Ordine','Q_LjungBox','GradiLiberta', ...
    'pValue','EsitoPortmanteau','RaggioSpettrale','IndiceStabilita', ...
    'EsitoStabilita','PercentConsistency'});

fprintf('\n=== CONFRONTO TRA ORDINE SOTTOSTIMATO E ORDINE CORRETTO ===\n');
disp(results);

%% 7. RAPPRESENTAZIONE GRAFICA
figure("Color", "w", "Name", "Validazione dell'ordine del VAR");
layout = tiledlayout(2,2, "TileSpacing", "compact", "Padding", "compact");

nexttile;
plot(1:h_test, residual_correlation_strength(1,:), "o-", ...
    "LineWidth", 1.4, "MarkerSize", 4);
hold on;
plot(1:h_test, residual_correlation_strength(2,:), "s-", ...
    "LineWidth", 1.4, "MarkerSize", 4);
grid on;
xlabel("Lag");
ylabel("||R_l||_F");
title("Correlazione residua complessiva");
legend("VAR(2): ordine imposto", ...
    sprintf("VAR(%d): ordine scelto con AIC", p_AIC), ...
    "Location", "best");

nexttile;
bar(model_orders, -log10(max(p_portmanteau, realmin)));
hold on;
yline(-log10(alpha), "r--", "Soglia \alpha = 0.05", ...
    "LineWidth", 1.2);
grid on;
xlabel("Ordine stimato");
ylabel("-log_{10}(p-value)");
title("Test portmanteau multivariato");

nexttile;
bar(model_orders, percent_consistency);
yline(85, "r--", "Soglia indicativa 85%", "LineWidth", 1.2);
grid on;
xlabel("Ordine stimato");
ylabel("Percent consistency [%]");
title("Consistenza del modello");

nexttile;
bar(model_orders, spectral_radius);
yline(1, "r--", "Limite di stabilita'", "LineWidth", 1.2);
grid on;
xlabel("Ordine stimato");
ylabel("Raggio spettrale");
title("Stabilita' del modello");

title(layout, ...
    "Effetto della scelta dell'ordine sulla validazione del modello VAR");

%% 8. MESSAGGIO DIDATTICO FINALE
fprintf('\nInterpretazione attesa:\n');
fprintf(['- Il VAR(2) omette le dipendenze ai ritardi 3-5. Il test portmanteau ' ...
         'dovrebbe quindi rifiutare H0.\n']);
fprintf(['- L''AIC seleziona p = %d senza conoscere a priori l''ordine reale. ' ...
         'I residui del modello selezionato dovrebbero risultare compatibili ' ...
         'con rumore bianco.\n'], p_AIC);
fprintf(['- La stabilita'' non garantisce da sola un buon adattamento: anche un ' ...
         'modello stabile puo'' avere ordine inadeguato e residui correlati.\n']);
```
