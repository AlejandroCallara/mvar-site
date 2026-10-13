%% HANDS-ON 2 - PDC, gPDC, DTF E dDTF IN UN SISTEMA VAR(2)
% Simulazione, stima ai minimi quadrati e analisi nel dominio della frequenza

% Convenzione degli indici:
%   misura(i,j,f) descrive la direzione j -> i alla frequenza f.

% Lo script calcola le versioni al quadrato di DTF, PDC, gPDC e dDTF:
%   DTF2_ij(f) = |H_ij(f)|^2 / sum_k |H_ik(f)|^2
%   PDC2_ij(f) = |A_ij(f)|^2 / sum_k |A_kj(f)|^2
%   gPDC2_ij(f): PDC pesata mediante le varianze delle innovazioni
%   dDTF2_ij(f) = ffDTF2_ij(f) * |P_ij(f)|^2

% Richiede MATLAB R2019b o successivo per tiledlayout.


clearvars;
close all;
clc;

rng(7, "twister");                 % Riproducibilita della simulazione

%% 1. PARAMETRI DELLA SIMULAZIONE
M = 4;                              % Numero di variabili
p = 2;                              % Ordine del modello VAR
N = 1000;                           % Numero di campioni mantenuti
dt = 1;                             % Periodo di campionamento [s]
fs = 1 / dt;                        % Frequenza di campionamento [Hz]
time = (0:N-1) * dt;

burnIn = 500;                       % Campioni eliminati come transitorio
rho = 0.9;                          % Modulo dei poli oscillatori
f1 = 0.1;                           % Frequenza propria di x1 [Hz]
f3 = 0.3;                           % Frequenza propria di x3 [Hz]

SigmaU = eye(M);                    % Covarianza delle innovazioni

%% 2. DEFINIZIONE DEL VAR(2) REALE
% x(t) = A1*x(t-1) + A2*x(t-2) + u(t)
Atrue = zeros(M, M, p);

% x1: oscillatore a 0.1 Hz
Atrue(1,1,1) = 2 * rho * cos(2*pi*f1);
Atrue(1,1,2) = -rho^2;

% Connessioni x1 -> x2 e x3 -> x2
Atrue(2,1,1) = 1;
Atrue(2,3,1) = 0.5;

% x3: oscillatore a 0.3 Hz e connessione x2 -> x3 a due ritardi
Atrue(3,3,1) = 2 * rho * cos(2*pi*f3);
Atrue(3,3,2) = -rho^2;
Atrue(3,2,1) = 0.5;
Atrue(3,2,2) = 0.5;

% Connessione x1 -> x4 al secondo ritardo
Atrue(4,1,2) = 1;

%% 3. GENERAZIONE DEI DATI SINTETICI
nTotal = N + burnIn + p;

% Queste realizzazioni standard e le condizioni iniziali saranno riutilizzate
% nella seconda simulazione, per isolare l'effetto di SigmaU.
standardInnovations = randn(M, nTotal);
initialState = randn(M, p);

U = chol(SigmaU, "lower") * standardInnovations;
X = zeros(M, nTotal);
X(:,1:p) = initialState;

for t = p+1:nTotal
    X(:,t) = Atrue(:,:,1) * X(:,t-1) ...
           + Atrue(:,:,2) * X(:,t-2) ...
           + U(:,t);
end

% Rimozione del transitorio iniziale
X = X(:, burnIn+p+1:end);

%% 4. STIMA MULTIVARIATA AI MINIMI QUADRATI
% Matrice delle risposte: ogni colonna contiene x(t)
Y = X(:, p+1:N);

% Matrice dei regressori impilati per ritardo:
% Z = [x(t-1); x(t-2)]
Z = zeros(M*p, N-p);
for lag = 1:p
    rows = (lag-1)*M + (1:M);
    Z(rows,:) = X(:, p+1-lag:N-lag);
end

% Y = B*Z + E, con B = [A1 A2].
% L'operatore / risolve il problema ai minimi quadrati senza calcolare
% esplicitamente inv(Z*Z').
Bhat = Y / Z;

Ahat = zeros(M, M, p);
for lag = 1:p
    columns = (lag-1)*M + (1:M);
    Ahat(:,:,lag) = Bhat(:,columns);
end

residuals = Y - Bhat * Z;
nEffective = size(Y,2);
SigmaUhat = residuals * residuals' / (nEffective - M*p);

fprintf("Raggio spettrale del VAR reale:  %.4f\n", ...
    spectralRadius(Atrue));
fprintf("Raggio spettrale del VAR stimato: %.4f\n\n", ...
    spectralRadius(Ahat));

disp("A1 reale:");
disp(Atrue(:,:,1));
disp("A1 stimata:");
disp(Ahat(:,:,1));
disp("A2 reale:");
disp(Atrue(:,:,2));
disp("A2 stimata:");
disp(Ahat(:,:,2));
disp("Covarianza stimata delle innovazioni:");
disp(SigmaUhat);

%% 5. RAPPRESENTAZIONE DEL VAR NEL DOMINIO DELLA FREQUENZA
% A(f) = I - sum_{r=1}^p A_r exp(-i*2*pi*f*r*dt)
% H(f) = A(f)^(-1)
nFrequencyBins = 512;
frequency = linspace(0, fs/2, nFrequencyBins);

Af = complex(zeros(M, M, nFrequencyBins));
Hf = complex(zeros(M, M, nFrequencyBins));

for n = 1:nFrequencyBins
    Af(:,:,n) = eye(M);

    for lag = 1:p
        phase = exp(-1i * 2*pi * frequency(n) * lag * dt);
        Af(:,:,n) = Af(:,:,n) - Ahat(:,:,lag) * phase;
    end

    % Numericamente preferibile a inv(Af(:,:,n))
    Hf(:,:,n) = Af(:,:,n) \ eye(M);
end

%% 6. CALCOLO DI DTF, PDC, gPDC E dDTF
DTF2 = zeros(M, M, nFrequencyBins);
PDC2 = zeros(M, M, nFrequencyBins);
gPDC2 = zeros(M, M, nFrequencyBins);
partialCoherence2 = zeros(M, M, nFrequencyBins);

% La gPDC usa le varianze stimate delle innovazioni.
innovationVariances = max(real(diag(SigmaUhat)), eps);

for n = 1:nFrequencyBins
    Hpower = abs(Hf(:,:,n)).^2;
    Apower = abs(Af(:,:,n)).^2;

    % Normalizzazione per riga: afflusso totale verso il bersaglio i
    DTF2(:,:,n) = Hpower ./ sum(Hpower, 2);

    % Normalizzazione per colonna: deflusso totale dalla sorgente j
    PDC2(:,:,n) = Apower ./ sum(Apower, 1);

    % gPDC: ciascuna riga i viene pesata mediante 1/SigmaUhat(i,i)
    weightedApower = Apower ./ innovationVariances;
    gPDC2(:,:,n) = weightedApower ./ sum(weightedApower, 1);

    % Matrice spettrale e coerenza parziale.
    % Se precisionSpectrum = S(f)^(-1), allora
    % P_ij(f) = -precisionSpectrum_ij / sqrt(P_ii * P_jj).
    spectralMatrix = Hf(:,:,n) * SigmaUhat * Hf(:,:,n)';
    precisionSpectrum = spectralMatrix \ eye(M);
    precisionDiagonal = max(real(diag(precisionSpectrum)), eps);
    normalization = sqrt(precisionDiagonal * precisionDiagonal.');
    partialCoherence = -precisionSpectrum ./ normalization;
    partialCoherence2(:,:,n) = abs(partialCoherence).^2;
end

% ffDTF: il denominatore contiene tutte le sorgenti e tutte le frequenze.
ffDTF2 = zeros(M, M, nFrequencyBins);
for target = 1:M
    fullFrequencyDenominator = sum(abs(Hf(target,:,:)).^2, "all");
    ffDTF2(target,:,:) = abs(Hf(target,:,:)).^2 ...
        / fullFrequencyDenominator;
end

dDTF2 = ffDTF2 .* partialCoherence2;

fprintf("Massima differenza PDC-gPDC con SigmaU = I: %.3e\n\n", ...
    max(abs(PDC2(:) - gPDC2(:))));

%% 7. RAPPRESENTAZIONE DEI SEGNALI E DELLE MISURE
samplesToShow = 1:min(200, N);

figure("Color", "w", "Name", "Dati sintetici");
signalLayout = tiledlayout(M, 1, "TileSpacing", "compact", ...
    "Padding", "compact");

for i = 1:M
    nexttile;
    plot(time(samplesToShow), X(i,samplesToShow), "LineWidth", 1);
    ylabel(sprintf("x_%d", i));
    grid on;
end
xlabel(signalLayout, "Tempo [s]");
title(signalLayout, "Prime 200 osservazioni del sistema VAR(2)");

plotConnectivityMeasure(frequency, DTF2, ...
    "DTF^2: flusso totale dalla sorgente j al bersaglio i");

plotConnectivityMeasure(frequency, PDC2, ...
    "PDC^2: connessioni dirette dalla sorgente j al bersaglio i");

% La scala viene scelta automaticamente perché la normalizzazione della
% ffDTF su tutte le frequenze produce valori molto inferiori a 1.
plotConnectivityMeasure(frequency, dDTF2, ...
    "dDTF^2: componente diretta del flusso direzionato", []);

%% 8. EFFETTO DELLE VARIANZE DELLE INNOVAZIONI: PDC vs gPDC
% Le dinamiche e i coefficienti del VAR non cambiano. Cambiano soltanto le
% varianze delle innovazioni, mantenendo le stesse realizzazioni standard.
SigmaUunequal = diag([1, 4, 0.25, 2]);
Uunequal = chol(SigmaUunequal, "lower") * standardInnovations;

Xunequal = zeros(M, nTotal);
Xunequal(:,1:p) = initialState;

for t = p+1:nTotal
    Xunequal(:,t) = Atrue(:,:,1) * Xunequal(:,t-1) ...
                  + Atrue(:,:,2) * Xunequal(:,t-2) ...
                  + Uunequal(:,t);
end

Xunequal = Xunequal(:, burnIn+p+1:end);

% Stima del secondo dataset
[AhatUnequal, SigmaUhatUnequal] = estimateVAR(Xunequal, p);
[AfUnequal, ~] = frequencyRepresentation(AhatUnequal, frequency, dt);
[PDC2unequal, gPDC2unequal] = calculatePDCandGPDC( ...
    AfUnequal, SigmaUhatUnequal);

disp("Covarianza teorica delle innovazioni nella seconda condizione:");
disp(SigmaUunequal);
disp("Covarianza stimata delle innovazioni nella seconda condizione:");
disp(SigmaUhatUnequal);

plotConnectivityComparison(frequency, PDC2unequal, gPDC2unequal, ...
    "Varianze delle innovazioni diverse: confronto PDC^2-gPDC^2");

%% FUNZIONI LOCALI
function radius = spectralRadius(A)
    [M, ~, p] = size(A);
    companion = [reshape(A, M, M*p); ...
                 eye(M*(p-1)), zeros(M*(p-1), M)];
    radius = max(abs(eig(companion)));
end

function plotConnectivityMeasure(frequency, measure, figureTitle, yLimits)
    M = size(measure, 1);

    if nargin < 4
        yLimits = [0, 1];
    elseif isempty(yLimits)
        maximumValue = max(measure(:));
        yLimits = [0, max(1.05 * maximumValue, eps)];
    end

    figure("Color", "w", "Name", figureTitle);
    layout = tiledlayout(M, M, "TileSpacing", "compact", ...
        "Padding", "compact");

    for target = 1:M
        for source = 1:M
            nexttile;
            plot(frequency, squeeze(measure(target,source,:)), ...
                "LineWidth", 1.2);
            xlim([frequency(1), frequency(end)]);
            ylim(yLimits);
            grid on;
            title(sprintf("%d \\rightarrow %d", source, target));

            if target == M
                xlabel("Frequenza [Hz]");
            end
            if source == 1
                ylabel("Intensita");
            end
        end
    end

    title(layout, figureTitle);
end

function [Ahat, SigmaUhat] = estimateVAR(X, p)
    [M, N] = size(X);
    Y = X(:,p+1:N);
    Z = zeros(M*p, N-p);

    for lag = 1:p
        rows = (lag-1)*M + (1:M);
        Z(rows,:) = X(:,p+1-lag:N-lag);
    end

    Bhat = Y / Z;
    Ahat = zeros(M, M, p);

    for lag = 1:p
        columns = (lag-1)*M + (1:M);
        Ahat(:,:,lag) = Bhat(:,columns);
    end

    residuals = Y - Bhat * Z;
    SigmaUhat = residuals * residuals' / (size(Y,2) - M*p);
end

function [Af, Hf] = frequencyRepresentation(A, frequency, dt)
    [M, ~, p] = size(A);
    nFrequencyBins = numel(frequency);
    Af = complex(zeros(M, M, nFrequencyBins));
    Hf = complex(zeros(M, M, nFrequencyBins));

    for n = 1:nFrequencyBins
        Af(:,:,n) = eye(M);

        for lag = 1:p
            phase = exp(-1i * 2*pi * frequency(n) * lag * dt);
            Af(:,:,n) = Af(:,:,n) - A(:,:,lag) * phase;
        end

        Hf(:,:,n) = Af(:,:,n) \ eye(M);
    end
end

function [PDC2, gPDC2] = calculatePDCandGPDC(Af, SigmaUhat)
    [M, ~, nFrequencyBins] = size(Af);
    PDC2 = zeros(M, M, nFrequencyBins);
    gPDC2 = zeros(M, M, nFrequencyBins);
    innovationVariances = max(real(diag(SigmaUhat)), eps);

    for n = 1:nFrequencyBins
        Apower = abs(Af(:,:,n)).^2;
        PDC2(:,:,n) = Apower ./ sum(Apower, 1);

        weightedApower = Apower ./ innovationVariances;
        gPDC2(:,:,n) = weightedApower ./ sum(weightedApower, 1);
    end
end

function plotConnectivityComparison(frequency, firstMeasure, ...
    secondMeasure, figureTitle)

    M = size(firstMeasure, 1);
    figure("Color", "w", "Name", figureTitle);
    layout = tiledlayout(M, M, "TileSpacing", "compact", ...
        "Padding", "compact");

    for target = 1:M
        for source = 1:M
            nexttile;
            plot(frequency, squeeze(firstMeasure(target,source,:)), ...
                "LineWidth", 1.2);
            hold on;
            plot(frequency, squeeze(secondMeasure(target,source,:)), ...
                "--", "LineWidth", 1.2);
            hold off;

            if target == 1 && source == 1
                legend("PDC^2", "gPDC^2", "Location", "best");
            end

            xlim([frequency(1), frequency(end)]);
            ylim([0, 1]);
            grid on;
            title(sprintf("%d \\rightarrow %d", source, target));

            if target == M
                xlabel("Frequenza [Hz]");
            end
            if source == 1
                ylabel("Intensita");
            end
        end
    end

    title(layout, figureTitle);
end
