%% DEMO 47 - SNR DIVERSO TRA I CANALI
% Un sistema realmente simmetrico puo apparire direzionale quando un canale
% contiene piu rumore di misura dell'altro.
%
% Convenzione: GC_j_to_i = log(varianza ridotta / varianza completa).
% Il codice usa solo funzioni MATLAB di base.

clearvars;
close all;
clc;
rng(7, "twister");

%% 1. VAR(2) SIMMETRICO DEL MATERIALE SUPPLEMENTARE
M = 2;
p = 2;
N = 50000;
burnIn = 500;

A1 = [0.55  0.025;
      0.025 0.55];
A2 = [-0.8 -0.1;
      -0.1 -0.8];

SigmaU = [1 0.3;
          0.3 1];

%% 2. GENERAZIONE DEL PROCESSO LATENTE
U = chol(SigmaU, "lower") * randn(M, N + burnIn + p);
Xlatent = zeros(M, N + burnIn + p);

for t = p+1:size(Xlatent,2)
    Xlatent(:,t) = A1 * Xlatent(:,t-1) + ...
                   A2 * Xlatent(:,t-2) + U(:,t);
end

Xlatent = Xlatent(:, burnIn+p+1:end);

%% 3. DUE CONDIZIONI DI MISURA
% Caso A: entrambi i canali sono osservati senza rumore additivo.
Xclean = Xlatent;

% Caso B: rumore di misura indipendente aggiunto soltanto al canale 1.
% La deviazione standard del rumore e uguale a quella del segnale 1:
% potenza del rumore = potenza del segnale, quindi SNR = 0 dB.
noiseStd = std(Xlatent(1,:));
Xnoisy = Xlatent;
Xnoisy(1,:) = Xnoisy(1,:) + noiseStd * randn(1,N);

%% 4. STIMA DELLA GC FORWARD E TIME-REVERSED
dataSets = {Xclean, Xnoisy};
conditionNames = {"SNR uguale", "Rumore sul canale 1"};
GCforward = zeros(2,2);          % colonne: 1->2, 2->1
GCreversed = zeros(2,2);

for condition = 1:2
    for directionOfTime = 1:2
        X = dataSets{condition};

        if directionOfTime == 2
            X = fliplr(X);        % inversione dell'ordine temporale
        end

        Y = X(:,p+1:end);
        Zfull = [X(:,p:end-1);
                 X(:,p-1:end-2)];

        Bfull = Y / Zfull;
        Efull = Y - Bfull * Zfull;

        % GC 1 -> 2: nel modello ridotto di x2 rimuovo il passato di x1.
        Zreduced = [X(2,p:end-1);
                    X(2,p-1:end-2)];
        bred = Y(2,:) / Zreduced;
        ered = Y(2,:) - bred * Zreduced;
        GC_1to2 = log(var(ered,1) / var(Efull(2,:),1));

        % GC 2 -> 1: nel modello ridotto di x1 rimuovo il passato di x2.
        Zreduced = [X(1,p:end-1);
                    X(1,p-1:end-2)];
        bred = Y(1,:) / Zreduced;
        ered = Y(1,:) - bred * Zreduced;
        GC_2to1 = log(var(ered,1) / var(Efull(1,:),1));

        if directionOfTime == 1
            GCforward(condition,:) = [GC_1to2 GC_2to1];
        else
            GCreversed(condition,:) = [GC_1to2 GC_2to1];
        end
    end
end

%% 5. ASIMMETRIA NETTA E TIME-REVERSED GC
% Asimmetria positiva: prevale 1 -> 2.
% Asimmetria negativa: prevale 2 -> 1.
netForward = GCforward(:,1) - GCforward(:,2);
netReversed = GCreversed(:,1) - GCreversed(:,2);

% Difference-based time-reversed GC (TRGC).
% Una direzione fisica ritardata cambia segno invertendo il tempo, quindi
% le due asimmetrie si sottraggono costruttivamente. Un imbalance dovuto al
% rumore di misura tende invece a conservarsi e viene attenuato.
TRGC = netForward - netReversed;

fprintf("\nGC forward\n");
fprintf("                         1 -> 2      2 -> 1\n");
fprintf("SNR uguale              %8.4f    %8.4f\n", ...
    GCforward(1,1), GCforward(1,2));
fprintf("Rumore sul canale 1     %8.4f    %8.4f\n", ...
    GCforward(2,1), GCforward(2,2));

fprintf("\nGC dopo inversione temporale\n");
fprintf("                         1 -> 2      2 -> 1\n");
fprintf("SNR uguale              %8.4f    %8.4f\n", ...
    GCreversed(1,1), GCreversed(1,2));
fprintf("Rumore sul canale 1     %8.4f    %8.4f\n", ...
    GCreversed(2,1), GCreversed(2,2));

fprintf("\nAsimmetria netta\n");
fprintf("                         forward     reversed       TRGC\n");
fprintf("SNR uguale              %8.4f     %8.4f     %8.4f\n", ...
    netForward(1), netReversed(1), TRGC(1));
fprintf("Rumore sul canale 1     %8.4f     %8.4f     %8.4f\n", ...
    netForward(2), netReversed(2), TRGC(2));

figure("Color","w","Name","Effetto di SNR diversi sulla GC");

subplot(1,3,1);
bar(GCforward);
set(gca,"XTickLabel",conditionNames);
ylabel("GC");
legend("1 -> 2","2 -> 1","Location","best");
title("Serie originali");
grid on;

subplot(1,3,2);
bar(GCreversed);
set(gca,"XTickLabel",conditionNames);
ylabel("GC");
legend("1 -> 2","2 -> 1","Location","best");
title("Serie invertite nel tempo");
grid on;

subplot(1,3,3);
bar([netForward netReversed TRGC]);
set(gca,"XTickLabel",conditionNames);
yline(0,"k-");
ylabel("Asimmetria direzionale");
legend("net forward","net reversed","TRGC","Location","best");
title("Correzione mediante time reversal");
grid on;

