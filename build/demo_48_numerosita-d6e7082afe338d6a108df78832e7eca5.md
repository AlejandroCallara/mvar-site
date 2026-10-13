# Demo 48 - Numerosità campionaria

## Scopo

Una simulazione Monte Carlo stima la GC usando un numero crescente di epoche. Serve a distinguere la struttura reale della rete dall'incertezza introdotta da campioni insufficienti.

## Risultati attesi

- con poche epoche le distribuzioni delle stime sono larghe;
- aumentando i dati diminuisce la deviazione standard delle stime;
- un test significativo non implica automaticamente assenza di bias.

I valori numerici cambiano con le realizzazioni casuali. Lo script imposta il generatore quando è necessario per rendere il confronto riproducibile.

## Output verificato

![Output MATLAB di Demo 48 - Numerosità campionaria](../assets/demo/demo_48_01.png)

## Esecuzione

Scaricare [Demo_48_numerosita_campionaria.m](../matlab/Demo_48_numerosita_campionaria.m), aprirlo in MATLAB ed eseguirlo per sezioni oppure integralmente. Gli script sono stati verificati con MATLAB R2025b.

## Codice MATLAB

```matlab
%% DEMO 48 - NUMEROSITA CAMPIONARIA E STIMA DELLA GC
% Lo stesso VAR(2) viene stimato usando un numero crescente di epoche.
% Con pochi dati la stima e piu variabile e puo essere distorta.
%
% Ogni epoca contiene 100 campioni, come nel materiale supplementare
% (durata 0.5 s, frequenza di campionamento 200 Hz).
% Il codice usa solo funzioni MATLAB di base.

clearvars;
close all;
clc;
rng(11, "twister");

%% 1. PARAMETRI
M = 2;
p = 2;
samplesPerTrial = 100;
burnIn = 200;
trialCounts = [5 10 50 100 500];
nMonteCarlo = 100;

A1 = [0.55  0.025;
      0.025 0.55];
A2 = [-0.8 -0.1;
      -0.1 -0.8];
SigmaU = [1 0.3;
          0.3 1];

GC_1to2 = zeros(nMonteCarlo, numel(trialCounts));
GC_2to1 = zeros(nMonteCarlo, numel(trialCounts));

%% 2. MONTE CARLO
for simulation = 1:nMonteCarlo
    for c = 1:numel(trialCounts)
        nTrials = trialCounts(c);

        % Y e Z raccolgono solo coppie temporali appartenenti alla stessa
        % epoca: non si creano regressori attraverso il confine tra epoche.
        Y = zeros(M, nTrials * (samplesPerTrial-p));
        Z = zeros(M*p, nTrials * (samplesPerTrial-p));
        firstColumn = 1;

        for trial = 1:nTrials
            U = chol(SigmaU,"lower") * ...
                randn(M, samplesPerTrial + burnIn + p);
            X = zeros(M, samplesPerTrial + burnIn + p);

            for t = p+1:size(X,2)
                X(:,t) = A1 * X(:,t-1) + ...
                         A2 * X(:,t-2) + U(:,t);
            end

            X = X(:,burnIn+p+1:end);
            columns = firstColumn:firstColumn + samplesPerTrial-p-1;
            Y(:,columns) = X(:,p+1:end);
            Z(:,columns) = [X(:,p:end-1); X(:,p-1:end-2)];
            firstColumn = columns(end) + 1;
        end

        % Modello completo bivariato.
        Bfull = Y / Z;
        Efull = Y - Bfull * Z;

        % Modello ridotto per 1 -> 2.
        Zreduced = Z([2 4],:);       % passato del solo canale 2
        bred = Y(2,:) / Zreduced;
        ered = Y(2,:) - bred * Zreduced;
        GC_1to2(simulation,c) = ...
            log(var(ered,1) / var(Efull(2,:),1));

        % Modello ridotto per 2 -> 1.
        Zreduced = Z([1 3],:);       % passato del solo canale 1
        bred = Y(1,:) / Zreduced;
        ered = Y(1,:) - bred * Zreduced;
        GC_2to1(simulation,c) = ...
            log(var(ered,1) / var(Efull(1,:),1));
    end
end

%% 3. MEDIA E DISPERSIONE DELLE STIME
mean12 = mean(GC_1to2,1);
mean21 = mean(GC_2to1,1);
std12 = std(GC_1to2,0,1);
std21 = std(GC_2to1,0,1);

fprintf("\nNumero di epoche, GC media e deviazione standard\n");
fprintf(" epoche       1 -> 2              2 -> 1\n");
for c = 1:numel(trialCounts)
    fprintf("%6d     %.4f +/- %.4f     %.4f +/- %.4f\n", ...
        trialCounts(c), mean12(c), std12(c), mean21(c), std21(c));
end

%% 4. GRAFICI
figure("Color","w","Name","Numerosita campionaria e GC");

subplot(1,2,1);
errorbar(trialCounts, mean12, std12, "o-", "LineWidth",1.3);
hold on;
errorbar(trialCounts, mean21, std21, "s-", "LineWidth",1.3);
set(gca,"XScale","log");
xlabel("Numero di epoche da 100 campioni");
ylabel("GC stimata");
legend("1 -> 2","2 -> 1","Location","best");
title("Media +/- deviazione standard");
grid on;

subplot(1,2,2);
hold on;
selectedCases = [1 3 5];
for k = 1:numel(selectedCases)
    c = selectedCases(k);
    horizontalJitter = 0.14 * (rand(nMonteCarlo,1) - 0.5);
    scatter(k + horizontalJitter, GC_1to2(:,c), 18, ...
        "filled", "MarkerFaceAlpha",0.35);
    plot([k-0.25 k+0.25], [mean12(c) mean12(c)], ...
        "k-", "LineWidth",2);
end
xlim([0.5 3.5]);
set(gca,"XTick",1:3,"XTickLabel",{"5","50","500"});
xlabel("Numero di epoche");
ylabel("GC 1 -> 2");
title("Singole stime; linea nera = media");
grid on;
```
