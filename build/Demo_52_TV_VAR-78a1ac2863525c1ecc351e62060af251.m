%% DEMO 52 - TV-VAR: finestre mobili, RLS, Kalman e PDC tempo-frequenza
% Un solo dataset x di dimensione M x N, senza trial.
% I tre stimatori osservano gli stessi campioni e usano lo stesso ordine.
%
% Rete simulata: 1 -> 2 -> 3, senza connessione diretta 1 -> 3.
% Variano nel tempo:
%   - il coefficiente a_21(1,t);
%   - il coefficiente a_32(1,t);
%   - la frequenza propria del primo oscillatore.
%
% Convenzione: A(i,j,r,t) e PDC(i,j,f,t) rappresentano j -> i.
% Si calcola la PDC QUADRATICA standard, normalizzata per colonna.
% Le mappe teoriche sono le PDC del modello "congelato" al tempo t.
% Non sono spettri globali di un processo stazionario.
%
% Tutto il codice e le funzioni locali sono contenuti in questo file.
% Richiede solo MATLAB base, versione R2022a o successiva (comando clim).
% Non cancella il workspace e non chiude eventuali figure gia' aperte.

%% 1. Configurazione riproducibile
rng(52,"twister");

fs = 200;                         % Hz
dt = 1/fs;
duration = 120;                    % secondi
N = round(duration*fs);
time = (0:N-1)/fs;
M = 3;
p = 2;                            % stesso ordine vero e stimato
d = M*p;                          % regressori per equazione
burnIn = 5*fs;

W = round(6*fs);                   % finestre causali di 6 secondi
step = round(0.10*fs);             % nuova mappa ogni 100 ms
initialSamples = round(8*fs);      % intervallo iniziale quasi stazionario
lambda = 1 - 2/W;                 % eta' media dei pesi circa W/2
q = 2e-7;                         % covarianza del random walk dei coefficienti
Qtheta = q*eye(M*d);

frequency = linspace(0,40,161);    % frequenze mostrate, in Hz
methodNames = ["Finestre","RLS","Kalman"];
methodColors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
evaluationSamples = initialSamples:step:N;
timeEvaluation = time(evaluationSamples);
nEvaluation = numel(evaluationSamples);

% Confronti eseguiti dopo l'assestamento e prima della fine del segnale.
validationMask = timeEvaluation >= 15 & timeEvaluation <= 115;

fprintf('\nDEMO 52: %d segnali, %d campioni, fs = %g Hz, VAR(%d)\n',M,N,fs,p);
fprintf('Finestre: W = %d, passo = %d campioni, assegnazione all''estremo destro.\n',W,step);
fprintf('RLS: lambda = %.6f, scala di memoria = %.2f s.\n',lambda,dt/(1-lambda));
fprintf('Kalman: Qtheta = %.2g I, Ru stimata dai primi %.1f s.\n',q,initialSamples/fs);

%% 2. Parametri veri: variazioni graduali e rete a catena
% Le transizioni sono smussate per non imporre una variazione istantanea.
rise21 = 0.5*(1+tanh((time-30)/3));
fall21 = 0.5*(1+tanh((time-78)/3));
coupling21 = 0.04 + 0.34*(rise21-fall21);
coupling32 = 0.34 - 0.28*0.5*(1+tanh((time-62)/4));
frequency1 = 10 + 5*0.5*(1+tanh((time-55)/6));
frequency2 = 22*ones(1,N);
frequency3 = 30*ones(1,N);

poleRadius = [0.90 0.85 0.80];
trueFrequency = [frequency1;frequency2;frequency3];
Atrue = zeros(M,M,p,N);
for channel = 1:M
    Atrue(channel,channel,1,:) = reshape( ...
        2*poleRadius(channel)*cos(2*pi*trueFrequency(channel,:)/fs),1,1,1,N);
    Atrue(channel,channel,2,:) = -poleRadius(channel)^2;
end
Atrue(2,1,1,:) = reshape(coupling21,1,1,1,N);
Atrue(3,2,1,:) = reshape(coupling32,1,1,1,N);

% Tutte le matrici sono triangolari inferiori. Le connessioni feedforward
% non spostano i poli dei tre oscillatori autonomi.
% Il primo oscillatore varia lentamente. La stabilita' "congelata" da sola
% non dimostra in generale la stabilita' di qualunque sistema tempo-variante.
Rtrue = eye(M);                    % innovazioni indipendenti e omoschedastiche
Ltrue = chol(Rtrue,"lower");

%% 3. Simulazione di una sola realizzazione M x tempo
% Il burn-in usa i coefficienti iniziali; poi parte il programma temporale.
nTotal = burnIn+N;
xAll = zeros(M,nTotal);
innovationsAll = Ltrue*randn(M,nTotal);
for sample = p+1:nTotal
    modelSample = max(1,sample-burnIn);
    for lag = 1:p
        xAll(:,sample) = xAll(:,sample) ...
            + Atrue(:,:,lag,modelSample)*xAll(:,sample-lag);
    end
    xAll(:,sample) = xAll(:,sample)+innovationsAll(:,sample);
end
x = xAll(:,burnIn+1:end);          % unico dataset da analizzare
assert(isequal(size(x),[M N]) && all(isfinite(x),"all"), ...
    'La simulazione deve produrre una matrice M x N finita.');

% Non si aggiunge rumore di misura: u_t e' l'innovazione del processo VAR.
% Nessuna standardizzazione tempo-dipendente: altererebbe i coefficienti
% veri con cui vogliamo confrontare le stime.

%% 4. Regressori e inizializzazione comune dai soli dati passati
Z = zeros(d,N);
for lag = 1:p
    rows = (lag-1)*M+(1:M);
    Z(rows,p+1:N) = x(:,p+1-lag:N-lag);
end
trainingSamples = p+1:initialSamples;
Z0 = Z(:,trainingSamples);
Y0 = x(:,trainingSamples);
B0 = Y0/Z0;                       % M x Mp, nessuna conoscenza dei parametri veri
residual0 = Y0-B0*Z0;
Ru = residual0*residual0'/(numel(trainingSamples)-d);
Ru = (Ru+Ru')/2;
G0inverse = (Z0*Z0')\eye(d);
assert(rcond(Ru)>1e-10 && rcond(Z0*Z0')>1e-10, ...
    'Inizializzazione mal condizionata.');

% Baseline descrittiva: un solo VAR costante, stimato sull'intero segnale.
% Serve a quantificare il vantaggio della descrizione tempo-variante.
% Non e' un metodo di previsione online.
Bconstant = x(:,p+1:N)/Z(:,p+1:N);

%% 5. Soluzione 1: finestre mobili causali
% In ciascuna finestra i coefficienti vengono stimati indipendentemente.
% I primi p campioni costruiscono i regressori: W-p osservazioni utilizzabili.
% Ogni regressore appartiene alla finestra; nessun campione futuro.
Bwindow = zeros(M,d,nEvaluation);
for index = 1:nEvaluation
    sample = evaluationSamples(index);
    observations = sample-W+p+1:sample;
    Bwindow(:,:,index) = x(:,observations)/Z(:,observations);
end

%% 6. Soluzione 2: RLS con dimenticanza esponenziale
% P e' Mp x Mp: inversa della matrice dei regressori pesata.
% Non e' la covarianza completa M^2p x M^2p del filtro di Kalman.
Brls = zeros(M,d,nEvaluation);
Brls(:,:,1) = B0;
Bcurrent = B0;
P = G0inverse;
nextIndex = 2;
predictionErrorRls = nan(M,N);
for sample = initialSamples+1:N
    z = Z(:,sample);
    errorPrior = x(:,sample)-Bcurrent*z;
    predictionErrorRls(:,sample) = errorPrior;
    gain = P*z/(lambda+z'*P*z);
    Bcurrent = Bcurrent+errorPrior*gain';
    P = (P-gain*z'*P)/lambda;
    P = (P+P')/2;                 % simmetrizzazione contro errori numerici
    if nextIndex <= nEvaluation && sample == evaluationSamples(nextIndex)
        Brls(:,:,nextIndex) = Bcurrent;
        nextIndex = nextIndex+1;
    end
end

%% 7. Soluzione 3: Kalman completo sullo stato dei coefficienti
% theta = vec([A1 ... Ap]), impilamento PER COLONNE.
% x_t = H_t theta_t + u_t, H_t = kron(z_t',I_M).
% theta_t = theta_(t-1) + w_t, Cov(w_t) = Qtheta.
% Ru resta fissata dopo l'inizializzazione: si isola il ruolo di Qtheta.
Bkalman = zeros(M,d,nEvaluation);
Bkalman(:,:,1) = B0;
theta = B0(:);
Ptheta = kron(G0inverse,Ru);
Istate = eye(M*d);
nextIndex = 2;
predictionErrorKalman = nan(M,N);
normalizedInnovation = nan(1,N);
for sample = initialSamples+1:N
    z = Z(:,sample);
    H = kron(z',eye(M));

    % Predizione, F = I.
    thetaPrior = theta;
    Pprior = Ptheta+Qtheta;

    % Innovazione del filtro: comprende anche l'incertezza sullo stato.
    nu = x(:,sample)-H*thetaPrior;
    S = H*Pprior*H'+Ru;
    K = (Pprior*H')/S;
    predictionErrorKalman(:,sample) = nu;
    normalizedInnovation(sample) = nu'*(S\nu);

    % Correzione. Forma di Joseph per la covarianza, numericamente robusta.
    theta = thetaPrior+K*nu;
    IKH = Istate-K*H;
    Ptheta = IKH*Pprior*IKH'+K*Ru*K';
    Ptheta = (Ptheta+Ptheta')/2;
    if nextIndex <= nEvaluation && sample == evaluationSamples(nextIndex)
        Bkalman(:,:,nextIndex) = reshape(theta,M,d);
        nextIndex = nextIndex+1;
    end
end

%% 8. PDC teorica e stimata sulla stessa griglia tempo-frequenza
% PDC(i,j,f,t) = |A_ij(f,t)|^2 / sum_k |A_kj(f,t)|^2.
% Normalizzazione per SORGENTE j (colonna). Non si usa la DTF.
% Il risultato e' PDC quadratica, compresa tra 0 e 1.
Btruth = reshape(Atrue(:,:,:,evaluationSamples),M,d,nEvaluation);
Bmethods = cat(4,Bwindow,Brls,Bkalman);
PDCtrue = localPdcSquared(Btruth,M,p,frequency,fs);
PDCestimated = zeros(M,M,numel(frequency),nEvaluation,3);
spectralRadiusTrue = localSpectralRadius(Btruth,M,p);
spectralRadiusEstimated = zeros(3,nEvaluation);
for method = 1:3
    PDCestimated(:,:,:,:,method) = localPdcSquared( ...
        Bmethods(:,:,:,method),M,p,frequency,fs);
    spectralRadiusEstimated(method,:) = localSpectralRadius( ...
        Bmethods(:,:,:,method),M,p);
end
PDCconstant = localPdcSquared(Bconstant,M,p,frequency,fs);

assert(max(spectralRadiusTrue)<1, ...
    'Il modello teorico locale deve essere stabile.');
assert(max(abs(sum(PDCtrue,1)-1),[],"all")<1e-10, ...
    'La PDC deve essere normalizzata per colonna.');

%% 9. Verifica quantitativa contro i parametri noti
% Confrontiamo soltanto i coefficienti fuori diagonale per la PDC di rete.
% Un singolo dataset permette una verifica della simulazione, non una
% caratterizzazione generale delle prestazioni degli stimatori.
offDiagonal = ~eye(M);
offDiagonal = repmat(offDiagonal,1,1,numel(frequency),sum(validationMask));
pdcTruthValidation = PDCtrue(:,:,:,validationMask);
trueCoefficientValidation = Btruth(:,:,validationMask);
coefficientRmse = zeros(3,1);
pdcRmse = zeros(3,1);
nullPdcMean = zeros(3,1);
stablePercent = zeros(3,1);
peakFrequencyMae = zeros(3,1);

% Periodi lontani dalle transizioni per verificare i livelli dei collegamenti.
activePlateau = timeEvaluation>=40 & timeEvaluation<=48;
inactivePlateau = timeEvaluation>=95 & timeEvaluation<=110;
g21active = zeros(3,1);
g21inactive = zeros(3,1);
g32early = zeros(3,1);
g32late = zeros(3,1);
expected21active = mean(coupling21(evaluationSamples(activePlateau)));
expected21inactive = mean(coupling21(evaluationSamples(inactivePlateau)));
expected32early = mean(coupling32(evaluationSamples(activePlateau)));
expected32late = mean(coupling32(evaluationSamples(inactivePlateau)));

constantError = Bconstant-trueCoefficientValidation;
constantCoefficientRmse = sqrt(mean(constantError.^2,"all"));
constantPdcError = PDCconstant-pdcTruthValidation;
constantPdcRmse = sqrt(mean(constantPdcError(offDiagonal).^2));

for method = 1:3
    coefficientError = Bmethods(:,:,validationMask,method)-trueCoefficientValidation;
    coefficientRmse(method) = sqrt(mean(coefficientError.^2,"all"));
    pdcError = PDCestimated(:,:,:,validationMask,method)-pdcTruthValidation;
    pdcRmse(method) = sqrt(mean(pdcError(offDiagonal).^2));
    nullPdcMean(method) = mean(PDCestimated(3,1,:,validationMask,method),"all");
    stablePercent(method) = 100*mean(spectralRadiusEstimated(method,validationMask)<1);

    map21 = squeeze(PDCestimated(2,1,:,validationMask,method));
    [~,maximumBins] = max(map21,[],1);
    [~,trueMaximumBins] = max(squeeze(PDCtrue(2,1,:,validationMask)),[],1);
    peakFrequencyMae(method) = mean(abs(frequency(maximumBins)-frequency(trueMaximumBins)));

    g21active(method) = mean(Bmethods(2,1,activePlateau,method),"all");
    g21inactive(method) = mean(Bmethods(2,1,inactivePlateau,method),"all");
    g32early(method) = mean(Bmethods(3,2,activePlateau,method),"all");
    g32late(method) = mean(Bmethods(3,2,inactivePlateau,method),"all");
end

summary = table(methodNames',coefficientRmse,pdcRmse,nullPdcMean, ...
    stablePercent,peakFrequencyMae, ...
    'VariableNames',{'Metodo','RMSE_coefficienti','RMSE_PDC_offdiag', ...
    'PDC_media_1verso3_assente','Percentuale_stabile','MAE_picco_Hz'});
plateauSummary = table(["Vero";methodNames'], ...
    [expected21active;g21active],[expected21inactive;g21inactive], ...
    [expected32early;g32early],[expected32late;g32late], ...
    'VariableNames',{'Metodo','a21_attivo','a21_tardivo','a32_pre','a32_post'});
disp(summary);
disp(plateauSummary);
fprintf('VAR costante: RMSE coefficienti = %.4f; RMSE PDC off-diag = %.4f.\n', ...
    constantCoefficientRmse,constantPdcRmse);

% Verifiche didattiche, fissate prima del confronto con le stime.
% Non sono test di significativita' della connettivita'.
assert(all(coefficientRmse<constantCoefficientRmse), ...
    'Uno stimatore non migliora la baseline costante sui coefficienti.');
assert(all(pdcRmse<constantPdcRmse), ...
    'Uno stimatore non migliora la baseline costante sulla PDC.');
assert(all(stablePercent>=99), ...
    'Troppe stime localmente instabili: riesaminare ordine e adattamento.');
assert(all(nullPdcMean<0.05), ...
    'PDC media eccessiva sulla connessione diretta assente 1 -> 3.');
assert(all(peakFrequencyMae<2), ...
    'Il picco della PDC 1 -> 2 non segue bene il riferimento teorico.');
assert(all(abs(g21active-expected21active)<0.08) && ...
    all(abs(g21inactive-expected21inactive)<0.08) && ...
    all(abs(g32early-expected32early)<0.08) && ...
    all(abs(g32late-expected32late)<0.08), ...
    'I livelli dei collegamenti nei plateau non riproducono la simulazione.');
fprintf('Verifiche della demo superate: coefficienti, PDC, stabilita'' e controllo nullo.\n');

%% 10. Figura 1: segnali e programma dei parametri veri
figure('Name','Demo 52 - Simulazione','Color','w','Position',[50 50 1200 800]);
tiledlayout(3,2,'TileSpacing','compact');
signalMask = time<=8;
for channel = 1:M
    nexttile(2*channel-1);
    plot(time(signalMask),x(channel,signalMask),'k');
    xlabel('Tempo (s)'); ylabel(sprintf('x_%d',channel));
    title(sprintf('Segnale %d: primi 8 s',channel));
    grid on;
end
nexttile(2); plot(time,coupling21,'k','LineWidth',1.8);
xlabel('Tempo (s)'); ylabel('a_{21}(1,t)'); title('Connessione 1 \rightarrow 2'); grid on;
nexttile(4); plot(time,coupling32,'k','LineWidth',1.8);
xlabel('Tempo (s)'); ylabel('a_{32}(1,t)'); title('Connessione 2 \rightarrow 3'); grid on;
nexttile(6); plot(time,frequency1,'k','LineWidth',1.8);
xlabel('Tempo (s)'); ylabel('Hz'); title('Frequenza propria di x_1'); grid on;
sgtitle('Una sola realizzazione: rete 1 \rightarrow 2 \rightarrow 3');

%% 11. Figura 2: coefficienti teorici e tre stime, senza riallineamenti
figure('Name','Demo 52 - Coefficienti','Color','w','Position',[80 80 1300 800]);
tiledlayout(2,2,'TileSpacing','compact');
coefficientPairs = [2 1;3 2;1 1;3 1];
coefficientTitles = ["1 \rightarrow 2","2 \rightarrow 3", ...
    "Memoria di x_1: frequenza variabile","1 \rightarrow 3: connessione diretta assente"];
for panel = 1:4
    nexttile; hold on;
    target = coefficientPairs(panel,1);
    source = coefficientPairs(panel,2);
    plot(timeEvaluation,squeeze(Btruth(target,source,:)), ...
        'k','LineWidth',2.2,'DisplayName','Vero');
    for method = 1:3
        plot(timeEvaluation,squeeze(Bmethods(target,source,:,method)), ...
            'Color',methodColors(method,:),'LineWidth',1.1, ...
            'DisplayName',methodNames(method));
    end
    xlabel('Tempo (s)'); ylabel(sprintf('a_{%d%d}(1,t)',target,source));
    title(coefficientTitles(panel)); grid on; xlim([15 115]);
    if panel==1
        legend('Location','best');
    end
end
sgtitle('Coefficienti al ritardo 1: variabilita'' e ritardo di inseguimento');

%% 12. Figura 3: tre direzioni, verita' e tre metodi sulla stessa scala
figure('Name','Demo 52 - PDC tempo-frequenza','Color','w','Position',[30 40 1500 950]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
directions = [2 1;3 2;3 1];        % bersaglio, sorgente
columnNames = ["Vero",methodNames];
for direction = 1:3
    target = directions(direction,1);
    source = directions(direction,2);
    for column = 1:4
        nexttile;
        if column==1
            map = squeeze(PDCtrue(target,source,:,:));
        else
            map = squeeze(PDCestimated(target,source,:,:,column-1));
            % Non presentare una stima instabile come connettivita' valida.
            map(:,spectralRadiusEstimated(column-1,:)>=1) = NaN;
        end
        imagesc(timeEvaluation,frequency,map);
        axis xy; clim([0 1]); xlim([15 115]); ylim([0 40]);
        title(sprintf('%s: %d \\rightarrow %d',columnNames(column),source,target));
        xlabel('Tempo (s)');
        if column==1
            ylabel('Frequenza (Hz)');
        end
        if direction==1
            hold on;
            plot(timeEvaluation,frequency1(evaluationSamples),'w--','LineWidth',1);
        end
        if column==4
            colorbar;
        end
    end
end
colormap(parula);
sgtitle('PDC quadratica: scala comune [0,1]; tratteggio = frequenza propria di x_1');

%% 13. Figura 4: riepilogo quantitativo e stabilita'
figure('Name','Demo 52 - Verifica quantitativa','Color','w','Position',[100 100 1200 720]);
tiledlayout(2,2,'TileSpacing','compact');
nexttile;
comparisonLabels = ["VAR costante",methodNames];
comparisonCategories = categorical(comparisonLabels,comparisonLabels,'Ordinal',true);
bar(comparisonCategories,[constantCoefficientRmse;coefficientRmse]);
ylabel('RMSE'); title('Errore su tutti i coefficienti'); grid on;
nexttile;
bar(comparisonCategories,[constantPdcRmse;pdcRmse]);
ylabel('RMSE'); title('Errore PDC sulle direzioni fuori diagonale'); grid on;
nexttile; hold on;
plot(timeEvaluation,spectralRadiusTrue,'k','LineWidth',2,'DisplayName','Vero');
for method = 1:3
    plot(timeEvaluation,spectralRadiusEstimated(method,:), ...
        'Color',methodColors(method,:),'DisplayName',methodNames(method));
end
yline(1,'k--','Limite di stabilita''','HandleVisibility','off');
xlim([15 115]); xlabel('Tempo (s)'); ylabel('Raggio spettrale');
title('Stabilita'' del modello congelato'); grid on; legend('Location','best');
nexttile;
bar(categorical(methodNames,methodNames,'Ordinal',true),nullPdcMean);
ylabel('PDC media'); ylim([0 0.05]); grid on;
title('Controllo: 1 \rightarrow 3 diretto assente');
sgtitle('Confronto con i parametri noti: nessuna soglia di significativita''');

%% 14. Punti da discutere durante la lezione
% - Le stime sono causali. I cambiamenti vengono seguiti con un ritardo che
%   dipende da W, lambda e Qtheta. Non si riallineano le curve a posteriori.
% - La frequenza variabile sposta la regione di PDC elevata per 1 -> 2.
% - La PDC 1 -> 3 teorica e' nulla: esiste un cammino indiretto, non un arco.
% - Le mappe stimate contengono valori non nulli anche per archi assenti.
%   Questi valori non dimostrano connettivita' significativa.
% - Modificare W, lambda o q permette di vedere il compromesso tra rumore
%   della stima e capacita' di seguire i cambiamenti.
% - Parametri scelti per questa simulazione non sono valori universali.
% - Ru e' stimata dai dati iniziali; Rtrue non viene usata dagli stimatori.
% - Fissare p al valore vero isola l'effetto del metodo di inseguimento.

%% Funzioni locali: stesso file, nessuna dipendenza esterna
function pdc = localPdcSquared(B,M,p,frequency,fs)
% B: M x Mp x tempo. Il primo indice di pdc e' il bersaglio.
    nTime = size(B,3);
    nFrequency = numel(frequency);
    pdc = zeros(M,M,nFrequency,nTime);
    phase = exp(-1i*2*pi*(1:p)'*frequency/fs);
    for index = 1:nTime
        A = reshape(B(:,:,index),M,M,p);
        for bin = 1:nFrequency
            Af = eye(M);
            for lag = 1:p
                Af = Af-A(:,:,lag)*phase(lag,bin);
            end
            power = abs(Af).^2;
            pdc(:,:,bin,index) = power./sum(power,1);
        end
    end
end

function radius = localSpectralRadius(B,M,p)
% Autovalori della companion per ogni modello locale.
    nTime = size(B,3);
    radius = zeros(1,nTime);
    for index = 1:nTime
        companion = [B(:,:,index);eye(M*(p-1)) zeros(M*(p-1),M)];
        radius(index) = max(abs(eig(companion)));
    end
end
