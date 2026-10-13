clear
clc

% simulazione u(t)
N = 100;
mu=[0 0];
sigma=[0.3 0; 0 0.3];
%potete provare a inserire un valore di cxy
%con |cxy|<sigmax*sigmay
u=mvnrnd(mu, sigma, N)';

x1(1)=0;%le condizioni iniziali possono essere diverse
x2(1)=0;%se avessi un modello di ordine 2 servirebbero altre due condizioni
X = [x1;x2];
A = [-0.7 -0.2; -0.2 -0.7];
p = 1;
for k=2:N
    X(:,k)=A*X(:,k-1)+u(:,k);
end

%da ora in poi consideriamo le due serie temporali e verifichiamo i valori degli indici di Granger
%dovremmo stimare due modelli marginali e uno completo


% =========================================================================
% 1. MODELLO COMPLETO (VAR)
% =========================================================================
Y = X(:, p+1:end);
[M, T] = size(X);
num_obs = T - p;

% Costruzione di Z per il modello completo
Z = zeros(M*p + 1, num_obs);
Z(1, :) = 1;
for k = 1:p
    row_idx = (1 + (k-1)*M + 1) : (1 + k*M);
    Z(row_idx, :) = X(:, (p - k + 1):(T - k));
end

B = Y * Z' / (Z * Z'); % Nota: / è preferibile numericamente rispetto a *inv()
u_hat = Y - B * Z;

% Varianze dei residui del modello completo (sulla diagonale della covariance)
sigma1_full = var(u_hat(1, :));
sigma2_full = var(u_hat(2, :));


% =========================================================================
% 2. MODELLI MARGINALI (AR UNIVARIATI)
% =========================================================================

% --- Modello Marginal per x1 ---
Y1 = X(1, p+1:end); % Variabile target: x1 al tempo t

% Regressore per x1: costante 1 + solo i lag di x1
Z1 = zeros(p + 1, num_obs);
Z1(1, :) = 1;
for k = 1:p
    Z1(1 + k, :) = X(1, (p - k + 1):(T - k));
end

B1 = Y1 * Z1' / (Z1 * Z1');
u1_hat = Y1 - B1 * Z1;
sigma1_restricted = var(u1_hat); % Varianza residua marginale di x1


% --- Modello Marginale per x2 ---
Y2 = X(2, p+1:end); % Variabile target: x2 al tempo t

% Regressore per x2: costante 1 + solo i lag di x2
Z2 = zeros(p + 1, num_obs);
Z2(1, :) = 1;
for k = 1:p
    Z2(1 + k, :) = X(2, (p - k + 1):(T - k));
end

B2 = Y2 * Z2' / (Z2 * Z2');
u2_hat = Y2 - B2 * Z2;
sigma2_restricted = var(u2_hat); % Varianza residua marginale di x2


% =========================================================================
% 3. CALCOLO DEGLI INDICI DI CAUSALITÀ DI GRANGER
% =========================================================================
% Causalità da x2 verso x1 (F_{2->1}): quanto x2 aiuta a predire x1
F_2_to_1 = log(sigma1_restricted / sigma1_full);

% Causalità da x1 verso x2 (F_{1->2}): quanto x1 aiuta a predire x2
F_1_to_2 = log(sigma2_restricted / sigma2_full);

% Stampa dei risultati
fprintf('Indice di Granger x2 -> x1: %.4f\n', F_2_to_1);
fprintf('Indice di Granger x1 -> x2: %.4f\n', F_1_to_2);

% =========================================================================
% 3. LIKELIHOOD RATIO TEST (LR Test / Chi-Quadro)
% =========================================================================

N_obs = num_obs; % Numero di osservazioni effettive (T - p)
q = p;           % Numero di parametri vincolati per ciascuna equazione (lag di x2 o x1)

% --- LR TEST: x2 -> x1 ---
% Statistica LR = N * log( sigma^2_rest / sigma^2_full )
LR_stat_2_to_1 = N_obs * log(sigma1_restricted / sigma1_full);
p_val_2_to_1   = 1 - chi2cdf(LR_stat_2_to_1, q); % p-value dalla distribuzione Chi2(q)

% --- LR TEST: x1 -> x2 ---
LR_stat_1_to_2 = N_obs * log(sigma2_restricted / sigma2_full);
p_val_1_to_2   = 1 - chi2cdf(LR_stat_1_to_2, q); % p-value dalla distribuzione Chi2(q)


% =========================================================================
% 4. STAMPA DEI RISULTATI
% =========================================================================
alpha = 0.05; % Livello di significatività

fprintf('\n=== LIKELIHOOD RATIO TEST PER CAUSALITÀ DI GRANGER (alpha = %.2f) ===\n', alpha);

fprintf('\n1) Ipotesi H0: x2 non causa Granger x1\n');
fprintf('   Indice di Granger: %.4f\n', log(sigma1_restricted / sigma1_full));
fprintf('   Statistica LR (Chi2): %.4f\n', LR_stat_2_to_1);
fprintf('   Gradi di libertà (q): %d\n', q);
fprintf('   p-value:              %.4e\n', p_val_2_to_1);
if p_val_2_to_1 < alpha
    fprintf('   Risultato: RIFIUTO H0 -> x2 CAUSA x1 alla Granger\n');
else
    fprintf('   Risultato: NON RIFIUTO H0 -> Nessuna causalità da x2 a x1\n');
end

fprintf('\n2) Ipotesi H0: x1 non causa Granger x2\n');
fprintf('   Indice di Granger: %.4f\n', log(sigma2_restricted / sigma2_full));
fprintf('   Statistica LR (Chi2): %.4f\n', LR_stat_1_to_2);
fprintf('   Gradi di libertà (q): %d\n', q);
fprintf('   p-value:              %.4e\n', p_val_1_to_2);
if p_val_1_to_2 < alpha
    fprintf('   Risultato: RIFIUTO H0 -> x1 CAUSA x2 alla Granger\n');
else
    fprintf('   Risultato: NON RIFIUTO H0 -> Nessuna causalità da x1 a x2\n');
end


% =========================================================================
% 5. VALIDAZIONE: TEST PORTMANTEAU MULTIVARIATO DI LJUNG-BOX
% =========================================================================
% Il test verifica congiuntamente che le matrici di autocovarianza dei
% residui siano nulle fino al ritardo h:
%       H0: C_1 = C_2 = ... = C_h = 0.
% A differenza di lbqtest applicato separatamente ai singoli canali, questo
% test considera anche le cross-correlazioni tra residui diversi.

n_res = size(u_hat, 2);
h_portmanteau = min(10, floor(n_res/5));

if h_portmanteau <= p
    error('Il lag del test portmanteau deve essere maggiore dell''ordine p.');
end

% Autocovarianza residua a lag zero e sua inversa generalizzata
C0 = (u_hat * u_hat') / n_res;
C0_inv = pinv(C0);

Q_LB = 0;
for lag = 1:h_portmanteau
    % C_l = E[u(t)u(t-l)']; il denominatore n_res mantiene la stessa
    % convenzione per tutte le autocovarianze campionarie.
    C_lag = (u_hat(:, lag+1:end) * u_hat(:, 1:end-lag)') / n_res;

    Q_LB = Q_LB + trace(C_lag' * C0_inv * C_lag * C0_inv) ...
                    / (n_res - lag);
end
Q_LB = n_res * (n_res + 2) * Q_LB;

df_portmanteau = M^2 * (h_portmanteau - p);
p_portmanteau = 1 - chi2cdf(Q_LB, df_portmanteau);

fprintf('\n=== TEST PORTMANTEAU MULTIVARIATO ===\n');
fprintf('Lag massimo h:            %d\n', h_portmanteau);
fprintf('Statistica Q di Ljung-Box: %.4f\n', Q_LB);
fprintf('Gradi di libertà:          %d\n', df_portmanteau);
fprintf('p-value:                   %.4e\n', p_portmanteau);

if p_portmanteau < alpha
    fprintf('Risultato: RIFIUTO H0 -> i residui non sono bianchi.\n');
else
    fprintf('Risultato: NON RIFIUTO H0 -> residui compatibili con rumore bianco.\n');
end


% =========================================================================
% 6. VALIDAZIONE: STABILITA' DEL MODELLO VAR
% =========================================================================
% Estrazione delle matrici A_k dalla matrice dei coefficienti stimati B.
% La prima colonna di B contiene l'intercetta.
A_hat = zeros(M, M, p);
for lag = 1:p
    columns = 2 + (lag-1)*M : 1 + lag*M;
    A_hat(:, :, lag) = B(:, columns);
end

% Matrice companion di dimensione (M*p) x (M*p).
if p == 1
    A_companion = A_hat(:, :, 1);
else
    A_top = reshape(A_hat, M, M*p);
    A_bottom = [eye(M*(p-1)), zeros(M*(p-1), M)];
    A_companion = [A_top; A_bottom];
end

lambda_companion = eig(A_companion);
rho_spectral = max(abs(lambda_companion));
stability_index = log(rho_spectral);
is_stable = rho_spectral < 1;

fprintf('\n=== STABILITA'' DEL MODELLO VAR ===\n');
fprintf('Massimo modulo degli autovalori: %.4f\n', rho_spectral);
fprintf('Indice di stabilità log(rho):    %.4f\n', stability_index);

if is_stable
    fprintf('Risultato: modello STABILE (tutti gli autovalori hanno modulo < 1).\n');
else
    fprintf('Risultato: modello NON STABILE.\n');
end


% =========================================================================
% 7. VALIDAZIONE: PERCENT CONSISTENCY
% =========================================================================
% Si genera una serie della stessa lunghezza a partire dal VAR stimato e si
% confrontano le auto- e cross-correlazioni dei dati reali e simulati.
% PC = 100 * (1 - ||R_sim - R_real||_2 / ||R_real||_2)

h_consistency = h_portmanteau;
percent_consistency = NaN;

if is_stable
    burn_in = 500;
    n_sim_total = T + burn_in;

    % Covarianza ML delle innovazioni stimate.
    Sigma_u_hat = (u_hat * u_hat') / n_res;
    Sigma_u_hat = (Sigma_u_hat + Sigma_u_hat') / 2;

    % Generazione robusta di innovazioni gaussiane anche in presenza di
    % piccoli autovalori negativi dovuti all'arrotondamento numerico.
    [V_sigma, D_sigma] = eig(Sigma_u_hat);
    eigenvalues_sigma = max(diag(D_sigma), 0);
    innovation_factor = V_sigma * diag(sqrt(eigenvalues_sigma));
    u_sim = innovation_factor * randn(M, n_sim_total);

    X_sim = zeros(M, n_sim_total);
    X_sim(:, 1:p) = X(:, 1:p);
    intercept_hat = B(:, 1);

    for t_sim = p+1:n_sim_total
        X_sim(:, t_sim) = intercept_hat + u_sim(:, t_sim);
        for lag = 1:p
            X_sim(:, t_sim) = X_sim(:, t_sim) ...
                + A_hat(:, :, lag) * X_sim(:, t_sim-lag);
        end
    end

    X_sim = X_sim(:, burn_in+1:burn_in+T);

    % Centratura dei dati reali e simulati.
    X_real_centered = X - mean(X, 2);
    X_sim_centered = X_sim - mean(X_sim, 2);

    C0_real = (X_real_centered * X_real_centered') / T;
    C0_sim = (X_sim_centered * X_sim_centered') / T;

    std_real = sqrt(max(diag(C0_real), eps));
    std_sim = sqrt(max(diag(C0_sim), eps));
    normalization_real = std_real * std_real';
    normalization_sim = std_sim * std_sim';

    % Ogni colonna contiene la matrice di correlazione vettorializzata a un
    % determinato lag, da zero fino a h_consistency.
    R_real = zeros(M^2, h_consistency + 1);
    R_sim = zeros(M^2, h_consistency + 1);

    for lag = 0:h_consistency
        C_real_lag = (X_real_centered(:, lag+1:end) ...
                    * X_real_centered(:, 1:end-lag)') / (T-lag);
        C_sim_lag = (X_sim_centered(:, lag+1:end) ...
                   * X_sim_centered(:, 1:end-lag)') / (T-lag);

        corr_real_lag = C_real_lag ./ normalization_real;
        corr_sim_lag = C_sim_lag ./ normalization_sim;

        R_real(:, lag+1) = corr_real_lag(:);
        R_sim(:, lag+1) = corr_sim_lag(:);
    end

    R_real_vector = R_real(:);
    R_sim_vector = R_sim(:);
    percent_consistency = 100 * ...
        (1 - norm(R_sim_vector - R_real_vector, 2) ...
        / norm(R_real_vector, 2));

    fprintf('\n=== CONSISTENZA DEL MODELLO ===\n');
    fprintf('Lag massimo considerato: %d\n', h_consistency);
    fprintf('Percent consistency:      %.2f %%\n', percent_consistency);

    if percent_consistency > 85
        fprintf('Risultato: il modello riproduce adeguatamente la struttura di correlazione.\n');
    else
        fprintf('Risultato: consistenza inferiore alla soglia indicativa dell''85%%.\n');
    end
else
    fprintf('\n=== CONSISTENZA DEL MODELLO ===\n');
    fprintf('Non calcolata: la simulazione da un VAR non stabile non è affidabile.\n');
end
