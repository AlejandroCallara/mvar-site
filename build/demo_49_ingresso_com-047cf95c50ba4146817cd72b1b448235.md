# Demo 49 - Ingresso comune non osservato

## Materiale originale

Questa demo riprende gli esempi e le cautele interpretative presentati in:

Bastos, A. M., & Schoffelen, J.-M. (2016). A tutorial review of functional connectivity analysis methods and their interpretational pitfalls. *Frontiers in Systems Neuroscience*, 9, 175.

## Scopo

La demo simula una sorgente latente che raggiunge due nodi con ritardi diversi. L'analisi bivariata produce una connessione apparente; il modello trivariato può condizionare sulla sorgente comune.

## Risultati attesi

- senza il nodo 3 compare GC apparente tra i nodi 1 e 2;
- includendo il nodo 3 la GC condizionata 1 verso 2 si riduce;
- la correzione funziona perché in questa simulazione il confondente è osservato nel modello completo.

I valori numerici cambiano con le realizzazioni casuali. Lo script imposta il generatore quando è necessario per rendere il confronto riproducibile.

## Output verificato

![Output MATLAB di Demo 49 - Ingresso comune non osservato](../assets/demo/demo_49_01.png)

## Esecuzione

Scaricare [Demo_49_ingresso_comune_non_osservato.m](../matlab/Demo_49_ingresso_comune_non_osservato.m), aprirlo in MATLAB ed eseguirlo per sezioni oppure integralmente. Gli script sono stati verificati con MATLAB R2025b.

## Codice MATLAB

```matlab
%% DEMO 49 - INGRESSO COMUNE NON OSSERVATO
% Il nodo 3 guida il nodo 1 al ritardo 1 e il nodo 2 al ritardo 2.
% Non esiste alcuna connessione diretta tra 1 e 2.
% Se il nodo 3 viene escluso, la GC bivariata puo suggerire 1 -> 2.
%
% Il codice usa solo funzioni MATLAB di base.

clearvars;
close all;
clc;
rng(21, "twister");

%% 1. VAR(2) DEL MATERIALE SUPPLEMENTARE
M = 3;
p = 2;
N = 50000;
burnIn = 500;

% x(t) = A1*x(t-1) + A2*x(t-2) + u(t)
A1 = [0.55 0    0.25;
      0    0.55 0;
      0    0    0.55];             % nodo 3 -> nodo 1, ritardo 1

A2 = [-0.8 0     0;
       0  -0.8  -0.1;
       0   0    -0.8];             % nodo 3 -> nodo 2, ritardo 2

SigmaU = eye(M);

%% 2. GENERAZIONE DEI DATI
U = chol(SigmaU,"lower") * randn(M, N + burnIn + p);
X = zeros(M, N + burnIn + p);

for t = p+1:size(X,2)
    X(:,t) = A1 * X(:,t-1) + ...
             A2 * X(:,t-2) + U(:,t);
end

X = X(:,burnIn+p+1:end);

%% 3. ANALISI BIVARIATA: IL NODO 3 NON E OSSERVATO
X12 = X(1:2,:);
Y12 = X12(:,p+1:end);
Z12 = [X12(:,p:end-1); X12(:,p-1:end-2)];
B12 = Y12 / Z12;
E12 = Y12 - B12 * Z12;

% GC 1 -> 2
Zred = Z12([2 4],:);
bred = Y12(2,:) / Zred;
ered = Y12(2,:) - bred * Zred;
GC_biv_1to2 = log(var(ered,1) / var(E12(2,:),1));

% GC 2 -> 1
Zred = Z12([1 3],:);
bred = Y12(1,:) / Zred;
ered = Y12(1,:) - bred * Zred;
GC_biv_2to1 = log(var(ered,1) / var(E12(1,:),1));

%% 4. ANALISI CONDIZIONATA: IL NODO 3 E INCLUSO
Y = X(:,p+1:end);
Z = [X(:,p:end-1); X(:,p-1:end-2)];
Bfull = Y / Z;
Efull = Y - Bfull * Z;

% GC 1 -> 2 | 3: elimino soltanto i ritardi della sorgente 1.
Zred = Z([2 3 5 6],:);
bred = Y(2,:) / Zred;
ered = Y(2,:) - bred * Zred;
GC_cond_1to2 = log(var(ered,1) / var(Efull(2,:),1));

% GC 2 -> 1 | 3: elimino soltanto i ritardi della sorgente 2.
Zred = Z([1 3 4 6],:);
bred = Y(1,:) / Zred;
ered = Y(1,:) - bred * Zred;
GC_cond_2to1 = log(var(ered,1) / var(Efull(1,:),1));

% Connessioni reali dal nodo comune.
Zred = Z([1 2 4 5],:);             % rimuovo il passato del nodo 3
bred = Y(1,:) / Zred;
ered = Y(1,:) - bred * Zred;
GC_3to1 = log(var(ered,1) / var(Efull(1,:),1));

bred = Y(2,:) / Zred;
ered = Y(2,:) - bred * Zred;
GC_3to2 = log(var(ered,1) / var(Efull(2,:),1));

%% 5. RISULTATI
fprintf("\nSenza osservare il nodo 3\n");
fprintf("GC 1 -> 2 = %.4f\n", GC_biv_1to2);
fprintf("GC 2 -> 1 = %.4f\n", GC_biv_2to1);

fprintf("\nCondizionando sul nodo 3\n");
fprintf("GC 1 -> 2 | 3 = %.4f\n", GC_cond_1to2);
fprintf("GC 2 -> 1 | 3 = %.4f\n", GC_cond_2to1);
fprintf("GC 3 -> 1 | 2 = %.4f\n", GC_3to1);
fprintf("GC 3 -> 2 | 1 = %.4f\n", GC_3to2);

figure("Color","w","Name","Ingresso comune non osservato");

subplot(1,2,1);
bar([GC_biv_1to2 GC_biv_2to1]);
set(gca,"XTickLabel",{"1 -> 2","2 -> 1"});
ylabel("GC");
title({"Nodo 3 escluso", "compare una direzione apparente"});
grid on;

subplot(1,2,2);
bar([GC_cond_1to2 GC_cond_2to1 GC_3to1 GC_3to2]);
set(gca,"XTickLabel",{"1 -> 2 | 3","2 -> 1 | 3", ...
                      "3 -> 1 | 2","3 -> 2 | 1"});
xtickangle(25);
ylabel("GC condizionata");
title({"Nodo 3 incluso", "emergono le connessioni reali"});
grid on;
```
