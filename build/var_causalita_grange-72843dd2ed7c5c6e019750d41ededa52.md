# VAR e causalità di Granger

## 1. Dal problema biomedico alla connettività direzionale

Molti segnali biomedici sono acquisiti simultaneamente da più canali. Un EEG multicanale, per esempio, descrive l'attività elettrica osservata in più sedi o ricostruita in più sorgenti corticali. La domanda non riguarda soltanto la presenza di oscillazioni comuni. Spesso interessa stabilire se il passato di un processo contiene informazione utile per prevederne un altro.

Questa domanda appartiene allo studio della **connettività effettiva**, cioè delle dipendenze direzionali tra processi. Va distinta dalla connettività funzionale, che descrive associazioni statistiche senza attribuire necessariamente una direzione. Una correlazione elevata tra due segnali non stabilisce quale dei due anticipi l'altro e può derivare da una terza sorgente comune.

```{figure} ../assets/figures/tipi_connettivita.png
:alt: Confronto tra connettività strutturale, funzionale ed effettiva
:width: 95%
:name: fig-tipi-connettivita

Tre livelli di descrizione della connettività: struttura anatomica, dipendenza statistica e influenza direzionale. Schema tratto dalla presentazione del corso.
```

La causalità di Granger formalizza una nozione precisa e limitata di causalità:

> Il processo $x_j$ causa $x_i$ nel senso di Granger se il passato di $x_j$ migliora la previsione di $x_i$ oltre l'informazione già contenuta nel passato di $x_i$ e delle altre variabili incluse nel modello.

I due ingredienti sono quindi la **precedenza temporale** e l'**informazione predittiva aggiuntiva**. Questa definizione non dimostra da sola un meccanismo fisico o biologico. La conclusione dipende dai segnali osservati, dal campionamento, dal preprocessing e dall'adeguatezza del modello.

```{admonition} Video di approfondimento
:class: tip
Per consolidare l'interpretazione predittiva della causalità di Granger si può vedere [*Granger causality (prediction)* di Mike X Cohen](https://www.youtube.com/watch?v=XqsSB_vpHLs). Il video è una risorsa integrativa e non sostituisce le verifiche sul modello discusse in questo capitolo.
```

## 2. Dall'autoregressione al modello vettoriale

### 2.1 Modello autoregressivo univariato

Un modello autoregressivo di ordine $p$ rappresenta il valore corrente di una serie mediante una combinazione lineare dei suoi valori passati:

$$
x(t)=c+\sum_{r=1}^{p}a_r x(t-r)+u(t),
$$

dove $u(t)$ è l'innovazione, cioè la parte di $x(t)$ che il passato non riesce a prevedere. L'ordine $p$ determina la memoria del modello. Con frequenza di campionamento $F_s$, la memoria esplicita è circa $p/F_s$ secondi.

### 2.2 Modello VAR multivariato

Raccogliamo $M$ segnali nel vettore

$$
\mathbf{x}(t)=
\begin{bmatrix}
x_1(t)&\cdots&x_M(t)
\end{bmatrix}^{\mathsf T}.
$$

Un VAR di ordine $p$ è

$$
\mathbf{x}(t)=\mathbf{c}+\sum_{r=1}^{p}\mathbf{A}_r\mathbf{x}(t-r)+\mathbf{u}(t),
$$

con $\mathbf{A}_r\in\mathbb{R}^{M\times M}$ e innovazioni

$$
\mathbf{u}(t)\sim\mathcal{N}(\mathbf{0},\boldsymbol{\Sigma}_u).
$$

L'elemento $A_{r,ij}$ quantifica il contributo di $x_j(t-r)$ alla previsione di $x_i(t)$. Gli elementi diagonali descrivono la dinamica autoregressiva di ciascun canale. Gli elementi fuori diagonale descrivono dipendenze ritardate tra canali.

Per un VAR senza intercetta il numero di coefficienti dinamici è $M^2p$. Aumentare il numero di variabili o l'ordine fa quindi crescere rapidamente il numero di parametri e la quantità di dati necessaria.

### 2.3 Esempio bivariato VAR(1)

Nel caso più semplice:

$$
\begin{aligned}
x_1(t)&=a_{11}x_1(t-1)+a_{12}x_2(t-1)+u_1(t),\\
x_2(t)&=a_{21}x_1(t-1)+a_{22}x_2(t-1)+u_2(t).
\end{aligned}
$$

Con la convenzione adottata, $a_{21}$ descrive $1\rightarrow2$ e $a_{12}$ descrive $2\rightarrow1$. Un coefficiente non nullo indica che il passato della sorgente entra nel modello del bersaglio. La causalità di Granger richiede però un confronto formale tra modelli, non la sola ispezione del coefficiente.

## 3. Assunzioni del modello VAR

### 3.1 Stazionarietà debole

Un processo è debolmente stazionario quando:

- la media non dipende da $t$;
- la varianza è finita e costante;
- la covarianza tra $\mathbf{x}(t)$ e $\mathbf{x}(t-h)$ dipende dal ritardo $h$, non dall'istante $t$.

Trend, cambiamenti di varianza e transizioni tra stati fisiologici violano questa ipotesi. Prima della stima possono servire detrending, rimozione della media, trasformazioni o una formulazione tempo-variante.

### 3.2 Stabilità

La stabilità è una proprietà del modello. Un VAR stabile dissipa l'effetto delle condizioni iniziali e delle perturbazioni. La matrice compagna è

$$
\mathbf{A}_c=
\begin{bmatrix}
\mathbf{A}_1&\mathbf{A}_2&\cdots&\mathbf{A}_p\\
\mathbf{I}&\mathbf{0}&\cdots&\mathbf{0}\\
&\ddots&\ddots&\vdots\\
\mathbf{0}&&\mathbf{I}&\mathbf{0}
\end{bmatrix}.
$$

Il modello è stabile se

$$
\rho(\mathbf{A}_c)=\max_k|\lambda_k(\mathbf{A}_c)|<1.
$$

In modo equivalente, le radici del polinomio autoregressivo inverso devono trovarsi fuori dal cerchio unitario. Sotto le usuali ipotesi sulle innovazioni, un VAR stabile ammette una soluzione causale stazionaria.

### 3.3 Innovazioni

Le innovazioni dovrebbero essere serialmente incorrelate e avere covarianza finita. La matrice $\boldsymbol{\Sigma}_u$ può avere elementi fuori diagonale: ciò rappresenta dipendenza istantanea tra innovazioni, non una direzione ritardata. L'ipotesi di gaussianità facilita la verosimiglianza e alcuni test, ma la stima ai minimi quadrati resta interpretabile anche in condizioni più generali.

## 4. Stima multivariata ai minimi quadrati

Per $t=p+1,\ldots,T$, definiamo

$$
\mathbf{z}(t)=
\begin{bmatrix}
1\\
\mathbf{x}(t-1)\\
\vdots\\
\mathbf{x}(t-p)
\end{bmatrix},
\qquad
\mathbf{B}=
\begin{bmatrix}
\mathbf{c}&\mathbf{A}_1&\cdots&\mathbf{A}_p
\end{bmatrix}.
$$

Impilando i campioni:

$$
\mathbf{Y}=
\begin{bmatrix}
\mathbf{x}(p+1)&\cdots&\mathbf{x}(T)
\end{bmatrix},
$$

$$
\mathbf{Z}=
\begin{bmatrix}
\mathbf{z}(p+1)&\cdots&\mathbf{z}(T)
\end{bmatrix}.
$$

Il problema diventa

$$
\mathbf{Y}=\mathbf{B}\mathbf{Z}+\mathbf{E}.
$$

La stima ai minimi quadrati minimizza l'energia dei residui:

$$
\widehat{\mathbf{B}}
=\arg\min_{\mathbf{B}}\|\mathbf{Y}-\mathbf{B}\mathbf{Z}\|_F^2
=\mathbf{Y}\mathbf{Z}^{\mathsf T}
(\mathbf{Z}\mathbf{Z}^{\mathsf T})^{-1}.
$$

Con la pseudoinversa:

$$
\widehat{\mathbf{B}}=\mathbf{Y}\mathbf{Z}^{\dagger}.
$$

In MATLAB è preferibile usare gli operatori di divisione matriciale invece di calcolare esplicitamente l'inversa:

```matlab
Bhat = Y / Z;
residuals = Y - Bhat * Z;
SigmaUhat = residuals * residuals' / size(residuals, 2);
```

I minimi quadrati costituiscono il riferimento per un VAR stazionario a coefficienti costanti. Con campioni moderati, regressori quasi collineari o modelli grandi, la stima può diventare instabile. Metodi stepwise, regolarizzati o lattice come Vieira-Morf possono migliorare alcune condizioni numeriche, ma richiedono scelte e verifiche proprie.

## 5. Scelta dell'ordine

Un ordine troppo basso lascia dipendenza temporale nei residui. Un ordine troppo alto aumenta la varianza delle stime e può adattare il rumore. I criteri di informazione bilanciano adattamento e complessità.

```{figure} ../assets/figures/ordine_modello.png
:alt: Compromesso nella scelta dell'ordine del modello autoregressivo
:width: 80%
:name: fig-ordine-modello

Compromesso nella scelta dell'ordine: un modello troppo semplice lascia struttura nei residui, mentre un modello troppo complesso aumenta varianza e rischio di sovradattamento. Schema tratto dalla presentazione del corso.
```

Indicando con $T_e=T-p$ il numero di campioni effettivi, con $k$ il numero di parametri e con $\widehat{\boldsymbol{\Sigma}}_u(p)$ la covarianza residua, una forma generale è

$$
IC(p)=\log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\text{penalità}(p,T_e).
$$

Forme comuni, a costanti additive irrilevanti, sono:

$$
\begin{aligned}
AIC(p)&=\log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{2k}{T_e},\\
BIC(p)&=\log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{k\log T_e}{T_e},\\
HQ(p)&=\log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{2k\log\log T_e}{T_e}.
\end{aligned}
$$

L'AIC applica una penalizzazione più moderata, il BIC tende a modelli più parsimoniosi, mentre HQ assume in genere un comportamento intermedio. Si sceglie il minimo del criterio nell'intervallo di ordini esaminato. Se la curva decresce fino al limite superiore, occorre ampliare l'intervallo o rivedere preprocessing e stazionarietà. Il criterio del gomito può fornire una guida euristica, non sostituisce la validazione.

L'ordine deve anche risultare fisiologicamente plausibile. Il ritardo atteso aiuta a definire l'intervallo di ricerca, ma non basta a identificare $p$.

## 6. Validazione del modello

Un valore di ordine selezionato da AIC o BIC non garantisce che il modello descriva adeguatamente i dati. La validazione combina più controlli.

### 6.1 Bianchezza dei residui

Se il VAR ha catturato la struttura temporale, le autocorrelazioni dei residui dovrebbero essere vicine a zero. Una prima ispezione usa le bande approssimate

$$
\pm\frac{2}{\sqrt{T_e}}.
$$

Circa il 95% dei coefficienti campionari dovrebbe cadere nelle bande sotto un'ipotesi di bianchezza. Questo controllo è punto per punto e può risultare fuorviante come test congiunto multivariato.

I test Portmanteau aggregano le autocorrelazioni fino al ritardo $h$. Una forma multivariata di Ljung-Box è

$$
Q_h=T_e(T_e+2)\sum_{\ell=1}^{h}
\frac{\operatorname{tr}\!\left(
\widehat{\mathbf{C}}_{\ell}^{\mathsf T}
\widehat{\mathbf{C}}_0^{-1}
\widehat{\mathbf{C}}_{\ell}
\widehat{\mathbf{C}}_0^{-1}
\right)}{T_e-\ell},
$$

dove $\widehat{\mathbf{C}}_{\ell}$ è l'autocovarianza residua al ritardo $\ell$. Una approssimazione comune usa

$$
Q_h\sim\chi^2_{M^2(h-p)}.
$$

Un $p$-value basso indica evidenza contro la bianchezza. Un $p$-value alto non dimostra che i residui siano bianchi.

### 6.2 Consistenza

La percent consistency confronta la struttura di correlazione dei dati reali con quella di dati simulati dal VAR stimato. Una forma è

$$
PC=100\left(1-
\frac{\|\mathbf{R}_{sim}-\mathbf{R}_{real}\|_2}
{\|\mathbf{R}_{real}\|_2}
\right).
$$

Valori elevati indicano che il modello riproduce bene la struttura correlativa complessiva. Non garantiscono la correttezza di ogni relazione direzionale.

### 6.3 Stabilità e controlli finali

Si verifica nuovamente $\rho(\widehat{\mathbf{A}}_c)<1$. Un modello può essere stabile ma inadeguato, per esempio quando l'ordine è troppo basso. Conviene inoltre controllare distribuzione e varianza dei residui, valori influenti, sensibilità all'ordine e robustezza a segmenti o prove diverse.

La [Demo 45](../demo/demo_45_validazione.md) mostra come un VAR stabile con ordine sottostimato possa fallire i test sui residui.

## 7. Causalità di Granger nel dominio del tempo

### 7.1 Modello completo e modello ridotto

Per testare $x_j\rightarrow x_i$, il **modello completo** include il passato di tutte le variabili. Il **modello ridotto** esclude tutti i ritardi della sorgente $x_j$ dall'equazione del bersaglio $x_i$.

Indichiamo con $\widehat{\sigma}_{i,F}^2$ e $\widehat{\sigma}_{i,R}^2$ le varianze residue del modello completo e ridotto. La causalità di Granger logaritmica è

$$
F_{j\rightarrow i}=\log\frac{\widehat{\sigma}_{i,R}^2}
{\widehat{\sigma}_{i,F}^2}.
$$

Poiché il modello ridotto ha meno regressori, ci aspettiamo $\widehat{\sigma}_{i,R}^2\geq\widehat{\sigma}_{i,F}^2$ e quindi $F_{j\rightarrow i}\geq0$ a meno di errori numerici o differenze nelle procedure di stima.

Nel caso multivariato la GC è **condizionata** sulle altre variabili incluse. Per testare $x_1\rightarrow x_2\mid x_3$, si rimuove il passato di $x_1$ dall'equazione di $x_2$, mantenendo quello di $x_2$ e $x_3$.

### 7.2 Dipendenza istantanea

Se $\Sigma_{u,ij}\neq0$, le innovazioni di $x_i$ e $x_j$ sono contemporaneamente correlate. Questa componente non è assegnabile a una direzione con i soli ritardi del VAR. Può riflettere interazioni più rapide del campionamento, ingressi comuni o mixing dei segnali.

Una misura simmetrica della dipendenza istantanea bivariata è

$$
F_{i\cdot j}=\log\frac{\Sigma_{u,ii}\Sigma_{u,jj}}
{\det\begin{bmatrix}
\Sigma_{u,ii}&\Sigma_{u,ij}\\
\Sigma_{u,ji}&\Sigma_{u,jj}
\end{bmatrix}}.
$$

### 7.3 Inferenza parametrica

I modelli completo e ridotto sono annidati. Un test $F$ può confrontare la riduzione della somma dei quadrati residui tenendo conto dei $p$ vincoli imposti dalla rimozione della sorgente. Sotto ipotesi regolari si possono anche usare test di Wald o likelihood-ratio.

Questi test richiedono residui e modello adeguati. Con campioni piccoli, non gaussianità marcata o procedure complesse di preprocessing, l'approssimazione asintotica può essere poco affidabile.

### 7.4 Bootstrap sotto l'ipotesi nulla

Una procedura empirica per $x_j\rightarrow x_i$ è:

1. stimare il modello ridotto che impone l'assenza della connessione;
2. ricampionare i vettori residui, mantenendo la loro struttura multivariata istantanea;
3. simulare nuovi dataset sotto $H_0$;
4. ristimare modello completo e ridotto per ogni dataset;
5. confrontare la GC osservata con la distribuzione nulla.

Con $B$ bootstrap, il $p$-value corretto è

$$
p=\frac{1+\#\{F_b\geq F_{obs}\}}{B+1}.
$$

La [Demo 50](../demo/demo_50_bootstrap.md) implementa questa procedura nelle due direzioni.

## 8. Limiti sperimentali e interpretativi

### 8.1 SNR diverso tra i canali

La GC confronta errori di predizione. Se un canale contiene più rumore di misura, la sua predicibilità cambia e l'asimmetria stimata può non riflettere l'accoppiamento reale. Anche un sistema simmetrico può apparire direzionale.

La time-reversed Granger causality ripete l'analisi sulle serie invertite nel tempo. Una dinamica fisica ritardata dovrebbe cambiare orientamento. Un'asimmetria dovuta soprattutto a rumore o mixing può invece mantenere lo stesso verso. La differenza tra asimmetria forward e reversed può attenuare questo effetto, ma non sostituisce un modello esplicito del rumore.

La [Demo 47](../demo/demo_47_snr.md) confronta due canali con SNR uguale e con rumore aggiunto a uno dei canali.

### 8.2 Numerosità campionaria

Il numero di parametri è dell'ordine di $M^2p$. Quando i campioni effettivi sono pochi rispetto ai parametri, aumentano varianza, bias e instabilità. La significatività statistica non elimina automaticamente il bias di stima.

Nelle analisi multitrial conta il numero complessivo di osservazioni utilizzabili e la loro indipendenza effettiva. Il confronto tra condizioni richiede numerosità comparabile o procedure che tengano conto dello sbilanciamento.

La [Demo 48](../demo/demo_48_numerosita.md) usa una simulazione Monte Carlo per mostrare la dispersione della GC al variare del numero di epoche.

### 8.3 Ingresso comune non osservato

Consideriamo una sorgente $x_3$ che raggiunge $x_1$ e $x_2$ con ritardi diversi. Se $x_3$ non viene misurata, $x_1$ può anticipare $x_2$ e il VAR bivariato può attribuire a $x_1\rightarrow x_2$ una connessione apparente.

La GC condizionata recupera la struttura corretta soltanto se il confondente è osservato e incluso con un modello adeguato. Non esiste una correzione puramente statistica che garantisca di rimuovere tutte le sorgenti latenti.

La [Demo 49](../demo/demo_49_ingresso_comune.md) confronta il modello bivariato con quello condizionato sul terzo nodo.

### 8.4 Mixing istantaneo e conduzione di volume

EEG e MEG misurano miscele lineari di sorgenti. La conduzione di volume produce dipendenze a ritardo zero e può alterare le stime direzionali, in particolare quando preprocessing, riferimento o ricostruzione delle sorgenti introducono fasi o ritardi. Le analisi dovrebbero dichiarare il livello dei segnali usati, il riferimento EEG e le procedure di leakage correction.

### 8.5 Preprocessing

Filtri non causali, downsampling, interpolazione di campioni e rimozione di artefatti possono modificare la struttura temporale. Le scelte di preprocessing devono essere identiche tra condizioni e documentate. È utile ripetere l'analisi con scelte ragionevoli alternative.

## 9. Dal VAR alla rappresentazione spettrale

Molti processi fisiologici mostrano oscillazioni. La rappresentazione in frequenza permette di stabilire in quali bande emerga una dipendenza direzionale.

Definiamo il polinomio autoregressivo

$$
\mathbf{A}(f)=\mathbf{I}-\sum_{r=1}^{p}
\mathbf{A}_r e^{-i2\pi f r\Delta t}.
$$

La trasformata del VAR soddisfa

$$
\mathbf{A}(f)\mathbf{X}(f)=\mathbf{U}(f).
$$

La funzione di trasferimento è

$$
\mathbf{H}(f)=\mathbf{A}(f)^{-1},
$$

da cui

$$
\mathbf{X}(f)=\mathbf{H}(f)\mathbf{U}(f).
$$

Se la covarianza delle innovazioni è $\boldsymbol{\Sigma}_u$, la matrice di densità spettrale è

$$
\mathbf{S}(f)=\mathbf{H}(f)\boldsymbol{\Sigma}_u\mathbf{H}^{\mathsf H}(f),
$$

dove $\mathsf H$ indica il trasposto coniugato. $S_{ii}(f)$ è lo spettro di $x_i$, mentre $S_{ij}(f)$ è il cross-spettro.

## 10. Coerenza e misure condizionate

La coerenza complessa è

$$
C_{ij}(f)=\frac{S_{ij}(f)}
{\sqrt{S_{ii}(f)S_{jj}(f)}},
\qquad 0\leq |C_{ij}(f)|^2\leq1.
$$

La coerenza misura associazione lineare in frequenza, ma non è direzionale. La fase può contenere informazione sui ritardi, però richiede cautela.

La coerenza immaginaria

$$
iCoh_{ij}(f)=\operatorname{Im}\{C_{ij}(f)\}
$$

riduce il contributo di accoppiamenti istantanei a fase nulla. Non elimina tutte le conseguenze del mixing.

Indicando con $\widetilde{\mathbf{S}}(f)=\mathbf{S}(f)^{-1}$ la precisione spettrale, la coerenza parziale è

$$
P_{ij}(f)=
-\frac{\widetilde{S}_{ij}(f)}
{\sqrt{\widetilde{S}_{ii}(f)\widetilde{S}_{jj}(f)}}.
$$

Questa misura descrive la coerenza residua tra $i$ e $j$ dopo il condizionamento lineare sulle altre variabili osservate.

## 11. Partial Directed Coherence

La PDC si costruisce dalla matrice $\mathbf{A}(f)$ e normalizza rispetto ai flussi uscenti dalla sorgente $j$:

$$
\pi_{i\leftarrow j}(f)=
\frac{A_{ij}(f)}
{\sqrt{\sum_{k=1}^{M}|A_{kj}(f)|^2}}.
$$

La versione quadratica soddisfa

$$
0\leq|\pi_{i\leftarrow j}(f)|^2\leq1,
\qquad
\sum_{i=1}^{M}|\pi_{i\leftarrow j}(f)|^2=1.
$$

La PDC è spesso interpretata come influenza diretta e condizionata da $j$ a $i$. La normalizzazione risponde alla domanda: come si distribuisce, tra i possibili bersagli, il flusso uscente dalla sorgente $j$?

```{figure} ../assets/figures/pdc_normalizzazione.png
:alt: Normalizzazione della PDC rispetto alle connessioni uscenti dalla sorgente
:width: 62%
:name: fig-pdc-normalizzazione

Nella PDC la connessione $j\rightarrow i$ evidenziata viene confrontata con tutte le connessioni uscenti dalla sorgente $j$. Schema tratto dalla presentazione del corso.
```

### 11.1 gPDC

La PDC generalizzata pesa le equazioni mediante la deviazione standard delle innovazioni:

$$
\overline{\pi}_{i\leftarrow j}(f)=
\frac{A_{ij}(f)/\sigma_i}
{\sqrt{\sum_{k=1}^{M}|A_{kj}(f)|^2/\sigma_k^2}},
\qquad \sigma_i^2=\Sigma_{u,ii}.
$$

Quando $\boldsymbol{\Sigma}_u=\mathbf{I}$, PDC e gPDC coincidono. Con scale o livelli di rumore diversi possono differire.

### 11.2 rPDC

La PDC rinormalizzata usa la covarianza della stima dei coefficienti VAR per produrre una statistica più adatta all'inferenza e meno dipendente dalle unità di misura. Il suo calcolo richiede la covarianza asintotica dei coefficienti e va distinto dalla sola normalizzazione geometrica della PDC.

## 12. Directed Transfer Function

La DTF usa la matrice di trasferimento e normalizza rispetto agli afflussi verso il bersaglio $i$:

$$
\gamma_{i\leftarrow j}(f)=
\frac{H_{ij}(f)}
{\sqrt{\sum_{k=1}^{M}|H_{ik}(f)|^2}}.
$$

La domanda è: quale quota dell'attività che raggiunge $i$ può essere associata alla sorgente $j$? Poiché $\mathbf{H}(f)$ contiene la propagazione attraverso l'intero sistema, la DTF può rilevare anche percorsi indiretti.

```{figure} ../assets/figures/dtf_normalizzazione.png
:alt: Normalizzazione della DTF rispetto alle connessioni entranti nel bersaglio
:width: 62%
:name: fig-dtf-normalizzazione

Nella DTF la connessione $j\rightarrow i$ evidenziata viene confrontata con tutti i contributi che raggiungono il bersaglio $i$. Schema tratto dalla presentazione del corso.
```

La full-frequency DTF usa un denominatore che integra sulle frequenze:

$$
\eta_{i\leftarrow j}^2(f)=
\frac{|H_{ij}(f)|^2}
{\sum_{\nu}\sum_{k=1}^{M}|H_{ik}(\nu)|^2}.
$$

La direct DTF combina ffDTF e coerenza parziale:

$$
\delta_{i\leftarrow j}^2(f)=
\eta_{i\leftarrow j}^2(f)|P_{ij}(f)|^2.
$$

Il fattore di coerenza parziale attenua i percorsi indiretti e mette in evidenza la componente diretta del flusso.

### 12.1 PDC e DTF rispondono a domande diverse

| Misura | Matrice di partenza | Normalizzazione | Interpretazione principale |
|---|---|---|---|
| PDC | $\mathbf{A}(f)$ | colonna, flussi uscenti dalla sorgente | connessione diretta e condizionata |
| gPDC | $\mathbf{A}(f)$ e $\boldsymbol{\Sigma}_u$ | uscita pesata per il rumore | PDC corretta per scale diverse |
| DTF | $\mathbf{H}(f)$ | riga, afflussi verso il bersaglio | trasferimento totale, anche indiretto |
| dDTF | $\mathbf{H}(f)$ e coerenza parziale | frequenze e condizionamento | componente diretta del trasferimento |

La scelta dipende dalla domanda. PDC e DTF non sono due stime intercambiabili della stessa quantità.

## 13. Inferenza per PDC e DTF

### 13.1 Distribuzioni asintotiche

La distribuzione dei coefficienti VAR può essere propagata alle misure spettrali con il metodo delta. Queste procedure risultano efficienti con campioni grandi e modello corretto, ma possono essere inaccurate vicino ai limiti dello spazio dei parametri, con campioni piccoli o con misure normalizzate non lineari.

### 13.2 Dati surrogati

I dati surrogati costruiscono empiricamente la distribuzione nulla. La procedura deve rispecchiare l'ipotesi che si vuole testare. Un surrogato che distrugge tutte le relazioni tra canali testa una null globale diversa dal bootstrap che annulla una sola connessione.

#### Phase randomization

Per ogni canale si calcola la FFT, si preserva il modulo e si randomizza la fase rispettando la simmetria coniugata. L'IFFT produce un segnale reale che conserva approssimativamente spettro di potenza e autocorrelazione del singolo canale, ma distrugge le relazioni di fase coerenti tra canali.

Per ciascun surrogato si ristima il VAR e si ricalcolano PDC e DTF. Il percentile $1-\alpha$ produce una soglia per direzione e frequenza. La [Demo 51](../demo/demo_51_surrogati_fase.md) mostra la procedura completa.

#### Causal shuffling

Il causal shuffling rompe selettivamente la relazione predittiva preservando altre caratteristiche dei segnali. La procedura esatta deve essere descritta perché differenti schemi di rimescolamento implicano differenti ipotesi nulle.

### 13.3 Confronti multipli

Una soglia puntuale al 95% non controlla l'errore complessivo quando si testano molte direzioni e frequenze. Opzioni comuni sono:

- controllo FDR sui $p$-value;
- max-statistic sui surrogati, usando il massimo sull'intera famiglia di test;
- test basati su cluster contigui in frequenza o tempo-frequenza;
- bande e finestre definite a priori.

Il numero di surrogati limita la risoluzione del $p$-value. Con $B$ surrogati, il minimo valore empirico con correzione $+1$ è $1/(B+1)$.

## 14. Esempio VAR(2) multivariato

Il secondo hands-on usa quattro processi:

$$
\begin{aligned}
x_1(t)&=2\rho\cos(2\pi f_1)x_1(t-1)-\rho^2x_1(t-2)+u_1(t),\\
x_2(t)&=x_1(t-1)+0.5x_3(t-1)+u_2(t),\\
x_3(t)&=2\rho\cos(2\pi f_3)x_3(t-1)-\rho^2x_3(t-2)\\
&\quad+0.5x_2(t-1)+0.5x_2(t-2)+u_3(t),\\
x_4(t)&=x_1(t-2)+u_4(t),
\end{aligned}
$$

con $\rho=0.9$, $f_1=0.1$ Hz, $f_3=0.3$ Hz e $\Delta t=1$ s. La struttura contiene $1\rightarrow2$, $3\rightarrow2$, $2\rightarrow3$ e $1\rightarrow4$. Il percorso $1\rightarrow2\rightarrow3$ può emergere nella DTF anche quando non esiste un coefficiente diretto $1\rightarrow3$.

```{figure} ../assets/figures/rete_var_faes_nollo.png
:alt: Rete VAR a quattro nodi dell'esempio di Faes e Nollo
:width: 70%
:name: fig-rete-var-faes-nollo

Struttura della rete VAR(2) usata nell'esempio. Adattata dal materiale di Faes e Nollo (2011) già incluso nella presentazione del corso.
```

Le due matrici di grafici seguenti mostrano perché PDC e DTF non devono essere lette come misure intercambiabili. Ogni pannello $(i,j)$ rappresenta la direzione $j\rightarrow i$.

```{figure} ../assets/figures/pdc_faes_nollo.png
:alt: PDC dell'esempio VAR a quattro nodi
:width: 90%
:name: fig-pdc-faes-nollo

PDC dell'esempio a quattro nodi: la normalizzazione privilegia le connessioni dirette uscenti dalla sorgente. Adattata da Faes e Nollo (2011).
```

```{figure} ../assets/figures/dtf_faes_nollo.png
:alt: DTF dell'esempio VAR a quattro nodi
:width: 90%
:name: fig-dtf-faes-nollo

DTF dello stesso esempio: la funzione di trasferimento può evidenziare anche la propagazione lungo percorsi indiretti. Adattata da Faes e Nollo (2011).
```

La [seconda esercitazione](../esercitazioni/esercitazione_2.md) guida stima e interpretazione.

## 15. Modelli VAR tempo-varianti

Un VAR stazionario usa coefficienti costanti. Nei dati evocati, nelle transizioni di stato o durante compiti, le interazioni possono cambiare:

$$
\mathbf{x}(t)=\sum_{r=1}^{p}\mathbf{A}_r(t)\mathbf{x}(t-r)+\mathbf{u}(t).
$$

### 15.1 Finestre scorrevoli

Si stima un VAR indipendente in ogni finestra. Il metodo è trasparente, ma assume stazionarietà locale. Una finestra lunga riduce la varianza e smussa le transizioni; una finestra corta migliora la localizzazione temporale ma rende la stima più instabile.

```{figure} ../assets/figures/lunghezza_finestra.png
:alt: Compromesso tra finestra temporale corta e lunga
:width: 82%
:name: fig-lunghezza-finestra

La lunghezza della finestra determina il compromesso tra risoluzione temporale e affidabilità della stima. Schema tratto dalla presentazione del corso.
```


### 15.2 Recursive Least Squares

L'RLS aggiorna i coefficienti campione per campione. In forma compatta:

$$
\mathbf{k}_t=
\frac{\mathbf{P}_{t-1}\mathbf{z}_t}
{\lambda+\mathbf{z}_t^{\mathsf T}\mathbf{P}_{t-1}\mathbf{z}_t},
$$

$$
\widehat{\mathbf{B}}_t=
\widehat{\mathbf{B}}_{t-1}+\mathbf{e}_t\mathbf{k}_t^{\mathsf T},
$$

$$
\mathbf{P}_t=\frac{1}{\lambda}
\left(\mathbf{P}_{t-1}-
\mathbf{k}_t\mathbf{z}_t^{\mathsf T}\mathbf{P}_{t-1}\right).
$$

Il fattore di dimenticanza $0<\lambda\leq1$ controlla il compromesso tra memoria e adattamento. In questo corso l'RLS serve soprattutto come termine di confronto; l'approfondimento riguarda il filtro di Kalman.

## 16. Filtro di Kalman per il TV-VAR

### 16.1 Modello stato-spazio generale

Il filtro di Kalman distingue uno stato nascosto e un'osservazione:

$$
\boldsymbol{\theta}_t=
\mathbf{F}_t\boldsymbol{\theta}_{t-1}+\mathbf{w}_t,
\qquad
\mathbf{w}_t\sim\mathcal{N}(\mathbf{0},\mathbf{Q}_{\theta}),
$$

$$
\mathbf{y}_t=
\mathbf{H}_t\boldsymbol{\theta}_t+\mathbf{v}_t,
\qquad
\mathbf{v}_t\sim\mathcal{N}(\mathbf{0},\mathbf{R}_y).
$$

Lo stato evolve secondo $\mathbf{F}_t$ e la matrice $\mathbf{Q}_{\theta}$ descrive quanto può cambiare. $\mathbf{H}_t$ collega lo stato alle osservazioni; $\mathbf{R}_y$ descrive il rumore osservativo.

Il filtro conserva sia la stima $\widehat{\boldsymbol{\theta}}$ sia la sua covarianza d'errore $\mathbf{P}$.

### 16.2 Predizione e correzione

**Predizione dello stato**

$$
\widehat{\boldsymbol{\theta}}_{t|t-1}=
\mathbf{F}_t\widehat{\boldsymbol{\theta}}_{t-1|t-1},
$$

$$
\mathbf{P}_{t|t-1}=
\mathbf{F}_t\mathbf{P}_{t-1|t-1}\mathbf{F}_t^{\mathsf T}
+\mathbf{Q}_{\theta}.
$$

**Innovazione**

$$
\boldsymbol{\nu}_t=
\mathbf{y}_t-\mathbf{H}_t\widehat{\boldsymbol{\theta}}_{t|t-1},
$$

$$
\mathbf{S}_t=
\mathbf{H}_t\mathbf{P}_{t|t-1}\mathbf{H}_t^{\mathsf T}
+\mathbf{R}_y.
$$

**Guadagno di Kalman**

$$
\mathbf{K}_t=
\mathbf{P}_{t|t-1}\mathbf{H}_t^{\mathsf T}\mathbf{S}_t^{-1}.
$$

**Correzione**

$$
\widehat{\boldsymbol{\theta}}_{t|t}=
\widehat{\boldsymbol{\theta}}_{t|t-1}
+\mathbf{K}_t\boldsymbol{\nu}_t,
$$

$$
\mathbf{P}_{t|t}=
(\mathbf{I}-\mathbf{K}_t\mathbf{H}_t)\mathbf{P}_{t|t-1}.
$$

Il guadagno bilancia la previsione e il nuovo campione. Se l'incertezza prevista è grande, il filtro corregge maggiormente lo stato. Se il rumore osservativo è grande, assegna meno peso al nuovo campione.

### 16.3 Il TV-VAR come modello stato-spazio

Lo stato raccoglie tutti i coefficienti:

$$
\boldsymbol{\theta}_t=
\operatorname{vec}
\begin{bmatrix}
\mathbf{A}_1(t)&\cdots&\mathbf{A}_p(t)
\end{bmatrix},
\qquad
\dim(\boldsymbol{\theta}_t)=M^2p.
$$

Definiamo

$$
\mathbf{z}_t=
\begin{bmatrix}
\mathbf{x}_{t-1}\\
\vdots\\
\mathbf{x}_{t-p}
\end{bmatrix}
\in\mathbb{R}^{Mp}.
$$

L'equazione di osservazione diventa

$$
\mathbf{x}_t=
\underbrace{(\mathbf{z}_t^{\mathsf T}\otimes\mathbf{I}_M)}_{\mathbf{H}_t}
\boldsymbol{\theta}_t+\mathbf{u}_t,
$$

con $\mathbf{H}_t\in\mathbb{R}^{M\times M^2p}$. Un modello comune per i coefficienti è il random walk:

$$
\boldsymbol{\theta}_t=
\boldsymbol{\theta}_{t-1}+\mathbf{w}_t,
$$

quindi $\mathbf{F}_t=\mathbf{I}$.

### 16.4 Ruolo delle covarianze

Se $\mathbf{Q}_{\theta}=q\mathbf{I}$, $q$ controlla la velocità di adattamento:

- $q$ piccolo produce coefficienti regolari ma può attenuare transizioni rapide;
- $q$ grande consente cambiamenti veloci ma può inseguire il rumore.

$\mathbf{R}_u=\operatorname{Cov}(\mathbf{u}_t)$ descrive l'incertezza dell'equazione VAR e può contenere dipendenze istantanee tra canali. Può essere fissata, inizializzata con un VAR stazionario o aggiornata ricorsivamente.

Nei dati multitrial, a ogni istante si possono usare più realizzazioni dello stesso processo, assumendo che condividano le stesse matrici $\mathbf{A}_r(t)$. La standardizzazione tra le prove deve essere descritta con precisione: media e deviazione standard vengono calcolate lungo la dimensione delle prove, separatamente per ogni tempo e canale, se questa è la procedura scelta.

### 16.5 Dalla stima alle mappe tempo-frequenza

A ogni istante si ricostruiscono $\widehat{\mathbf{A}}_r(t)$ e quindi

$$
\widehat{\mathbf{A}}(f,t)=
\mathbf{I}-\sum_{r=1}^{p}
\widehat{\mathbf{A}}_r(t)e^{-i2\pi fr\Delta t},
$$

$$
\widehat{\mathbf{H}}(f,t)=
\widehat{\mathbf{A}}(f,t)^{-1}.
$$

Si ottengono GC tempo-variante, PDC e DTF tempo-frequenza. La risoluzione formale può essere campione per campione, ma quella effettiva dipende da $\mathbf{Q}_{\theta}$, $\mathbf{R}_u$, inizializzazione, regolarizzazione e quantità di dati.

## 17. Inferenza sulle mappe tempo-frequenza

Una mappa contiene molti test correlati. L'interpretazione di singoli pixel sopra una soglia puntuale è fragile. Quando appropriato si usano max-statistic o test basati su cluster, mantenendo la struttura spazio-temporale nei surrogati.

Prima dell'analisi conviene definire bande e finestre temporali motivate dalla domanda fisiologica. Vanno riportati ordine del VAR, metodo adattivo, parametri di regolarizzazione, numero di surrogati e famiglia di confronti controllata.

## 18. Checklist per un'analisi riproducibile

Riportare almeno:

1. segnali, sorgenti, frequenza di campionamento e preprocessing;
2. riferimento EEG, numero di prove e campioni effettivi;
3. ordine del modello, metodo di stima e criteri di selezione;
4. controlli di stabilità e bianchezza dei residui;
5. definizione esatta e normalizzazione di GC, PDC o DTF;
6. modello nullo, numero di bootstrap o surrogati e correzione multipla;
7. analisi di sensibilità a ordine, finestra, parametri adattivi e preprocessing;
8. limiti legati a mixing, rumore e variabili non osservate.

```{admonition} Interpretazione finale
:class: warning
La causalità di Granger misura informazione predittiva direzionale condizionata sul modello e sulle variabili osservate. Una conclusione fisiologica richiede coerenza con il disegno sperimentale, robustezza analitica e controllo delle spiegazioni alternative.
```


