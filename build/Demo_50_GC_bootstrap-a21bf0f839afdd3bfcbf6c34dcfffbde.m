%% DEMO 1 - SIGNIFICATIVITA DELLA GC CON SURROGATI
% Tutti i passaggi sono svolti direttamente nello script.
%
% Convenzione: A(i,j) rappresenta la direzione j -> i.

clearvars;
close all;
clc;

rng(7);

%% Parametri
N = 1000;
burnIn = 500;
B = 500;                 % Aumentare a 1000 per l'analisi finale
alpha = 0.05;

% Sistema con connessione reale 1 -> 2
A = [-0.7,  0.0; ...
     -0.2, -0.7];

Sigma = [0.3, 0.0; ...
         0.0, 0.3];

%% Simulazione del VAR(1)
nTotal = N + burnIn;
u = chol(Sigma,"lower") * randn(2,nTotal);
x = zeros(2,nTotal);
x(:,1) = randn(2,1);

for t = 2:nTotal
    x(:,t) = A*x(:,t-1) + u(:,t);
end

% Eliminazione del transitorio
x = x(:,burnIn+1:end);

figure("Color","w");
subplot(2,1,1);
plot(x(1,:));
ylabel("x_1");
grid on;

subplot(2,1,2);
plot(x(2,:));
ylabel("x_2");
xlabel("Campione");
grid on;

%% GC osservata nelle due direzioni
% direzione 1: 1 -> 2, quindi target=2 e source=1
% direzione 2: 2 -> 1, quindi target=1 e source=2
targets = [2,1];
sources = [1,2];
labels = ["1 -> 2","2 -> 1"];

GCobserved = zeros(1,2);

for direction = 1:2
    target = targets(direction);
    source = sources(direction);

    y = x(target,2:end);
    Zfull = x(:,1:end-1);

    % Modello completo
    betaFull = y / Zfull;
    residualFull = y - betaFull*Zfull;
    varianceFull = mean(residualFull.^2);

    % Modello ridotto: viene eliminata la sorgente
    remainingVariable = setdiff(1:2,source);
    Zreduced = x(remainingVariable,1:end-1);
    betaReduced = y / Zreduced;
    residualReduced = y - betaReduced*Zreduced;
    varianceReduced = mean(residualReduced.^2);

    GCobserved(direction) = log(varianceReduced/varianceFull);
end

fprintf("GC osservata 1 -> 2: %.4f\n",GCobserved(1));
fprintf("GC osservata 2 -> 1: %.4f\n",GCobserved(2));

%% Bootstrap dei residui sotto H0
GCnull = zeros(B,2);

for direction = 1:2
    target = targets(direction);
    source = sources(direction);

    Y = x(:,2:end);
    Z = x(:,1:end-1);

    % Stima del VAR completo
    Acomplete = Y/Z;

    % Modello nullo: il coefficiente source -> target viene imposto a zero
    A0 = Acomplete;
    remainingVariable = setdiff(1:2,source);
    beta0 = Y(target,:) / Z(remainingVariable,:);
    A0(target,:) = 0;
    A0(target,remainingVariable) = beta0;

    % Residui del modello nullo
    residuals0 = Y - A0*Z;
    residuals0 = residuals0 - mean(residuals0,2);
    numberOfResiduals = size(residuals0,2);

    for b = 1:B
        % Ricampionamento con reinserimento dei vettori residui
        sampledColumns = randi(numberOfResiduals,1,N+burnIn-1);
        uStar = residuals0(:,sampledColumns);

        % Generazione del dataset surrogato sotto H0
        xStar = zeros(2,N+burnIn);
        xStar(:,1) = x(:,randi(N));

        for t = 2:N+burnIn
            xStar(:,t) = A0*xStar(:,t-1) + uStar(:,t-1);
        end

        xStar = xStar(:,burnIn+1:end);

        % GC del dataset surrogato
        yStar = xStar(target,2:end);
        ZfullStar = xStar(:,1:end-1);

        betaFullStar = yStar/ZfullStar;
        errorFullStar = yStar - betaFullStar*ZfullStar;

        ZreducedStar = xStar(remainingVariable,1:end-1);
        betaReducedStar = yStar/ZreducedStar;
        errorReducedStar = yStar - betaReducedStar*ZreducedStar;

        GCnull(b,direction) = log(mean(errorReducedStar.^2) ...
            / mean(errorFullStar.^2));
    end
end

%% Soglie e p-value empirici
pValue = zeros(1,2);
threshold = zeros(1,2);

for direction = 1:2
    pValue(direction) = ...
        (1 + sum(GCnull(:,direction) >= GCobserved(direction))) / (B+1);

    sortedGC = sort(GCnull(:,direction));
    thresholdIndex = ceil((1-alpha)*B);
    threshold(direction) = sortedGC(thresholdIndex);
end

fprintf("\n1 -> 2: soglia = %.4f, p = %.4f\n", ...
    threshold(1),pValue(1));
fprintf("2 -> 1: soglia = %.4f, p = %.4f\n", ...
    threshold(2),pValue(2));

%% Grafici delle distribuzioni nulle
figure("Color","w");

for direction = 1:2
    subplot(1,2,direction);
    histogram(GCnull(:,direction),30,"Normalization","pdf");
    hold on;
    xline(GCobserved(direction),"r","LineWidth",2);
    xline(threshold(direction),"k--","LineWidth",1.5);
    hold off;
    grid on;
    xlabel("GC");
    ylabel("Densita");
    title(sprintf("%s, p = %.4f",labels(direction),pValue(direction)));
    legend("H_0","GC osservata","Soglia 95%");
end
