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

I tre livelli della figura rispondono a domande diverse. La connettività **strutturale** riguarda l'esistenza di un supporto anatomico; la connettività **funzionale** descrive il fatto che due attività varino insieme; la connettività **effettiva** cerca invece di caratterizzare un'influenza orientata. Passare da un livello all'altro non è automatico: due regioni anatomicamente connesse possono non mostrare un'interazione misurabile in una certa condizione, mentre due segnali funzionalmente associati possono esserlo a causa di una sorgente comune. La causalità di Granger fornisce una definizione operativa dell'ultimo livello, costruita sui dati osservati e non su una descrizione biofisica completa della rete (Bressler & Seth, 2011; Bastos & Schoffelen, 2016).

La causalità di Granger formalizza una nozione precisa e limitata di causalità:

> Il processo $x_j$ causa $x_i$ nel senso di Granger se il passato di $x_j$ migliora la previsione di $x_i$ oltre l'informazione già contenuta nel passato di $x_i$ e delle altre variabili incluse nel modello.

I due ingredienti sono quindi la **precedenza temporale** e l'**informazione predittiva aggiuntiva**. Questa definizione non dimostra da sola un meccanismo fisico o biologico. La conclusione dipende dai segnali osservati, dal campionamento, dal preprocessing e dall'adeguatezza del modello.

È importante soffermarsi sull'espressione “predittiva aggiuntiva”. Non basta osservare che un evento in $x_j$ precede un evento simile in $x_i$: occorre dimostrare che il passato di $x_j$ riduce l'incertezza sulla previsione di $x_i$ dopo aver già utilizzato tutta l'informazione disponibile nel modello. In un sistema bivariato questo significa condizionare sul passato del bersaglio; in un sistema multivariato significa condizionare anche sulle altre variabili osservate. La GC è quindi sempre relativa a un insieme informativo dichiarato.

La direzionalità non implica necessariamente un'unica freccia. In un circuito con feedback possono essere significative entrambe le direzioni, ma con intensità o contenuto spettrale differenti. Inoltre, la relazione stimata riguarda i segnali e la scala temporale analizzata: dire “$x_j$ causa $x_i$ nel senso di Granger” è più preciso che affermare semplicemente “$x_j$ causa $x_i$”.

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

L'idea è simile a una previsione a un passo: al tempo $t$ si dispone dei campioni fino a $t-1$ e si costruisce la migliore previsione lineare del campione successivo. I coefficienti $a_r$ descrivono quanto ciascun ritardo contribuisce alla previsione, mentre $u(t)$ raccoglie ciò che rimane imprevedibile. Se il modello è adeguato, l'innovazione non contiene più una struttura temporale sfruttabile; se invece nei residui resta autocorrelazione, una parte della dinamica non è stata modellata.

Un modello AR non deve essere letto come una semplice interpolazione della forma d'onda. È un modello generativo: una volta specificati coefficienti e distribuzione delle innovazioni, può produrre nuove realizzazioni con proprietà statistiche simili a quelle osservate. Questa interpretazione rende naturali sia la simulazione sia la validazione del modello.

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

Una lettura pratica delle matrici consiste nel fissare prima la riga, cioè il bersaglio. La riga $i$ raccoglie tutti i contributi necessari a prevedere $x_i$; la colonna $j$ raccoglie invece gli effetti ritardati che partono da $x_j$. Questa convenzione evita uno degli errori più frequenti nell'interpretazione di PDC e DTF, dove l'elemento $(i,j)$ continua a rappresentare la direzione $j\rightarrow i$.

Il passaggio dal caso univariato al VAR introduce un vantaggio essenziale: tutte le equazioni vengono stimate come parti dello stesso sistema. In questo modo una relazione tra due canali può essere valutata tenendo conto degli altri segnali inclusi. Il prezzo è un rapido aumento dei parametri e, quindi, della quantità di dati necessaria per ottenere stime affidabili (Hytti et al., 2006).

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

Per esempio, se $a_{21}\neq0$, conoscere $x_1(t-1)$ aiuta potenzialmente a prevedere $x_2(t)$. La parola “potenzialmente” è necessaria perché il coefficiente stimato può essere diverso da zero per variabilità campionaria, collinearità o specificazione inadeguata del modello. La GC valuta l'effetto congiunto di tutti i $p$ coefficienti associati alla sorgente e domanda se la loro rimozione peggiori in modo apprezzabile la previsione.

## 3. Assunzioni del modello VAR

### 3.1 Stazionarietà debole

Un processo è debolmente stazionario quando:

- la media non dipende da $t$;
- la varianza è finita e costante;
- la covarianza tra $\mathbf{x}(t)$ e $\mathbf{x}(t-h)$ dipende dal ritardo $h$, non dall'istante $t$.

Trend, cambiamenti di varianza e transizioni tra stati fisiologici violano questa ipotesi. Prima della stima possono servire detrending, rimozione della media, trasformazioni o una formulazione tempo-variante.

La stazionarietà va riferita all'intervallo di analisi. Un segnale EEG può essere non stazionario sull'intera registrazione e risultare approssimativamente stazionario in una finestra breve o in una fase sperimentale omogenea. Ridurre la finestra, però, diminuisce il numero di osservazioni disponibili: si guadagna localizzazione temporale ma si perde precisione. La scelta è quindi un compromesso modellistico, non un semplice passaggio di preprocessing.

È inoltre utile distinguere la stazionarietà dei dati dalla stabilità del modello stimato. Un VAR stabile genera, sotto le ipotesi usuali, un processo stazionario; ma trovare coefficienti stabili non dimostra che il segmento reale analizzato fosse omogeneo. Per questo la verifica delle radici deve essere affiancata all'ispezione dei segnali, delle statistiche locali e dei residui.

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

Le innovazioni non sono necessariamente “rumore strumentale”. Rappresentano tutto ciò che non è prevedibile con i ritardi e le variabili inserite: attività endogena, ingressi esterni non misurati, componenti non lineari e rumore di misura. Una covarianza residua fuori diagonale segnala che due equazioni ricevono variazioni contemporanee comuni, ma un VAR strettamente causale non può attribuire a tali variazioni una direzione all'interno dello stesso campione (Faes & Nollo, 2011).

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

Dal punto di vista geometrico, la stima cerca la combinazione dei regressori ritardati che proietta ogni campione corrente nello spazio generato dal passato. Tutte le equazioni condividono la stessa matrice dei regressori $\mathbf{Z}$, ma hanno bersagli diversi nelle righe di $\mathbf{Y}$. La stima congiunta preserva così una rappresentazione coerente della rete.

Con dati multitrial, i campioni validi delle diverse prove possono contribuire alla stessa regressione, purché le prove siano considerate realizzazioni comparabili dello stesso processo. Non si devono però creare ritardi artificiali concatenando la fine di una prova con l'inizio della successiva: i regressori vanno costruiti separatamente in ogni trial e solo dopo impilati. Questa distinzione diventa decisiva nelle finestre brevi, dove l'uso di più prove fornisce le osservazioni necessarie senza allungare la finestra temporale.

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

La selezione dell'ordine non è quindi una gara tra criteri numerici. Un flusso di lavoro ragionevole consiste nel: definire un intervallo compatibile con la frequenza di campionamento; calcolare più criteri; esaminare la zona in cui il miglioramento si stabilizza; stimare i modelli candidati; controllare stabilità e residui. Se AIC propone un ordine elevato ma i residui sono già bianchi con un modello più parsimonioso, la scelta più semplice può essere preferibile. Al contrario, un minimo del BIC non è accettabile se lascia una chiara struttura predittiva nei residui.

L'ordine ha anche un significato spettrale. Un modello troppo corto non riesce a rappresentare oscillazioni strette o ritardi distribuiti; uno troppo lungo può introdurre picchi instabili e connessioni spurie. Per questo Hytti et al. (2006) sottolineano che l'ordine ottimale per la previsione a un passo non è necessariamente quello che restituisce le funzioni di trasferimento più interpretabili.

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

Questi controlli rispondono a domande complementari. La stabilità chiede se la dinamica stimata sia matematicamente ammissibile; la bianchezza chiede se resti informazione temporale non spiegata; la consistenza chiede se il modello riproduca la struttura globale dei dati. Nessun indicatore, preso isolatamente, certifica il modello. La validazione è piuttosto una ricerca di incoerenze: se un controllo fallisce, l'interpretazione causale deve essere sospesa e il modello rivisto.

## 7. Causalità di Granger nel dominio del tempo

### 7.1 Modello completo e modello ridotto

Per testare $x_j\rightarrow x_i$, il **modello completo** include il passato di tutte le variabili. Il **modello ridotto** esclude tutti i ritardi della sorgente $x_j$ dall'equazione del bersaglio $x_i$.

Indichiamo con $\widehat{\sigma}_{i,F}^2$ e $\widehat{\sigma}_{i,R}^2$ le varianze residue del modello completo e ridotto. La causalità di Granger logaritmica è

$$
F_{j\rightarrow i}=\log\frac{\widehat{\sigma}_{i,R}^2}
{\widehat{\sigma}_{i,F}^2}.
$$

Poiché il modello ridotto ha meno regressori, ci aspettiamo $\widehat{\sigma}_{i,R}^2\geq\widehat{\sigma}_{i,F}^2$ e quindi $F_{j\rightarrow i}\geq0$ a meno di errori numerici o differenze nelle procedure di stima.

Il modello completo e quello ridotto formulano due domande controfattuali ben delimitate. Il completo chiede quanto bene si possa prevedere il bersaglio usando tutto il passato disponibile; il ridotto chiede cosa accada quando si rende indisponibile il passato della sorgente, lasciando invariato il resto dell'insieme informativo. Il logaritmo del rapporto tra le varianze rende la misura nulla quando le due previsioni sono equivalenti e crescente quando l'esclusione della sorgente produce un peggioramento.

La GC non misura direttamente l'ampiezza di una risposta evocata né il valore di un singolo coefficiente. Due connessioni con coefficienti numericamente simili possono produrre valori diversi di GC, perché la misura dipende anche dalla dinamica autoregressiva, dalla covarianza delle innovazioni e dalla ridondanza informativa con gli altri segnali.

Nel caso multivariato la GC è **condizionata** sulle altre variabili incluse. Per testare $x_1\rightarrow x_2\mid x_3$, si rimuove il passato di $x_1$ dall'equazione di $x_2$, mantenendo quello di $x_2$ e $x_3$.

Il condizionamento è il principale vantaggio dell'analisi multivariata, ma va interpretato con attenzione. Se $x_3$ media davvero il percorso $x_1\rightarrow x_3\rightarrow x_2$, condizionare su $x_3$ può eliminare l'effetto indiretto e isolare la componente diretta. Se invece $x_3$ è una misura rumorosa o una variabile collidente, la sua inclusione può modificare la stima in modo meno intuitivo. L'insieme delle variabili deve quindi essere motivato dalla fisiologia e dal disegno sperimentale, non scelto soltanto in base alla disponibilità dei canali.

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

Significatività ed entità dell'effetto sono due informazioni diverse. Con molti campioni anche una riduzione minima dell'errore può risultare statisticamente significativa; con pochi campioni un effetto rilevante può non superare la soglia. È quindi opportuno riportare insieme valore di GC, intervallo o distribuzione empirica, numerosità effettiva e procedura di correzione per confronti multipli.

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

Il punto concettuale è che la GC opera sugli errori di predizione dei segnali misurati, non direttamente sugli stati neurali latenti. Se il rumore di misura rende una serie meno prevedibile, cambia sia il denominatore sia il numeratore del rapporto di varianze. Bastos e Schoffelen (2016) mostrano così che un accoppiamento simmetrico può produrre una direzione prevalente apparente quando i canali hanno SNR diversi. Il controllo delle potenze e del rumore non è quindi un'analisi accessoria, ma parte dell'interpretazione della direzionalità.

### 8.2 Numerosità campionaria

Il numero di parametri è dell'ordine di $M^2p$. Quando i campioni effettivi sono pochi rispetto ai parametri, aumentano varianza, bias e instabilità. La significatività statistica non elimina automaticamente il bias di stima.

Nelle analisi multitrial conta il numero complessivo di osservazioni utilizzabili e la loro indipendenza effettiva. Il confronto tra condizioni richiede numerosità comparabile o procedure che tengano conto dello sbilanciamento.

La [Demo 48](../demo/demo_48_numerosita.md) usa una simulazione Monte Carlo per mostrare la dispersione della GC al variare del numero di epoche.

La numerosità introduce due problemi distinti. Il primo è la varianza delle stime VAR: con pochi dati, i coefficienti oscillano maggiormente tra realizzazioni. Il secondo è il bias delle misure normalizzate e non negative: anche sotto l'ipotesi nulla, una stima finita non si concentra esattamente a zero. Per confrontare soggetti o condizioni è quindi preferibile mantenere comparabili numero di trial e campioni, oppure costruire procedure di ricampionamento che rendano esplicito lo sbilanciamento.

### 8.3 Ingresso comune non osservato

Consideriamo una sorgente $x_3$ che raggiunge $x_1$ e $x_2$ con ritardi diversi. Se $x_3$ non viene misurata, $x_1$ può anticipare $x_2$ e il VAR bivariato può attribuire a $x_1\rightarrow x_2$ una connessione apparente.

La GC condizionata recupera la struttura corretta soltanto se il confondente è osservato e incluso con un modello adeguato. Non esiste una correzione puramente statistica che garantisca di rimuovere tutte le sorgenti latenti.

La [Demo 49](../demo/demo_49_ingresso_comune.md) confronta il modello bivariato con quello condizionato sul terzo nodo.

Il ritardo diverso con cui un ingresso comune raggiunge due nodi è particolarmente insidioso: crea una precedenza temporale reale nei segnali osservati, pur in assenza di una connessione diretta tra essi. Un modello bivariato non possiede l'informazione necessaria per distinguere questa geometria da una vera influenza $x_1\rightarrow x_2$. L'esempio chiarisce un limite generale: il condizionamento protegge dalle cause comuni soltanto quando esse sono misurate con qualità sufficiente e incluse nel modello.

### 8.4 Mixing istantaneo e conduzione di volume

EEG e MEG misurano miscele lineari di sorgenti. La conduzione di volume produce dipendenze a ritardo zero e può alterare le stime direzionali, in particolare quando preprocessing, riferimento o ricostruzione delle sorgenti introducono fasi o ritardi. Le analisi dovrebbero dichiarare il livello dei segnali usati, il riferimento EEG e le procedure di leakage correction.

Una dipendenza a fase zero non deve essere automaticamente scartata come artefatto, né considerata automaticamente fisiologica. Il mixing può distribuire la stessa sorgente su più sensori, ma interazioni biologiche rapide o ingressi comuni possono anch'essi produrre covariazioni quasi istantanee. Metodi che attenuano le componenti a fase zero riducono alcuni falsi positivi al prezzo di perdere una parte potenzialmente informativa del segnale. La scelta va quindi collegata alla domanda sperimentale e verificata con analisi di sensibilità.

### 8.5 Preprocessing

Filtri non causali, downsampling, interpolazione di campioni e rimozione di artefatti possono modificare la struttura temporale. Le scelte di preprocessing devono essere identiche tra condizioni e documentate. È utile ripetere l'analisi con scelte ragionevoli alternative.

In particolare, restringere una banda con un filtro e poi applicare una misura direzionale non equivale sempre a calcolare la causalità spettrale del modello broadband. Il filtro modifica autocorrelazioni, fasi e ordine richiesto dal VAR; filtri bidirezionali possono inoltre introdurre dipendenze temporali non compatibili con una lettura causale ingenua. Quando l'obiettivo è localizzare la GC in frequenza, è spesso preferibile stimare un modello adeguato sui dati opportunamente campionati e decomporre successivamente la misura nel dominio spettrale (Barnett & Seth, 2014).

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

La trasformazione in frequenza non introduce una nuova dinamica: riorganizza la stessa informazione contenuta nei coefficienti temporali. $\mathbf{A}(f)$ descrive, frequenza per frequenza, i vincoli diretti imposti dal modello; $\mathbf{H}(f)$ descrive come le innovazioni si propagano attraverso l'intera rete. Questa dualità spiega perché PDC e DTF, pur partendo dallo stesso VAR, enfatizzino aspetti differenti.

La rappresentazione spettrale è particolarmente utile nei segnali fisiologici perché una relazione può essere debole se integrata su tutte le frequenze ma concentrata in una banda. Un picco direzionale va tuttavia letto insieme agli spettri di potenza: una misura normalizzata può assumere valori elevati anche dove la potenza assoluta è piccola, e quindi dove il contributo fisiologico complessivo è limitato.

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

Il termine “diretta” va inteso rispetto alle variabili contenute nel modello. La PDC è nulla per $j\rightarrow i$ quando tutti i coefficienti VAR fuori diagonale corrispondenti sono nulli; un percorso mediato da un nodo osservato non crea da solo una PDC diretta. Se il mediatore non viene osservato, però, il modello può rappresentare il percorso latente mediante un legame apparente tra i nodi disponibili.

La normalizzazione per colonna rende ogni sorgente una sorta di bilancio relativo: se una sorgente distribuisce la propria influenza tra molti bersagli, la quota assegnata a ciascuno può ridursi anche senza una diminuzione del coefficiente assoluto. Perciò la PDC è adatta a confrontare la distribuzione delle uscite di una sorgente, ma i suoi valori non vanno interpretati come quantità assolute di potenza trasferita (Faes & Nollo, 2011).

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

La DTF segue il segnale lungo la rete. Se $1\rightarrow2$ e $2\rightarrow3$, un'innovazione generata nel nodo 1 può contribuire allo spettro del nodo 3 anche quando non esiste alcun coefficiente diretto $1\rightarrow3$. Questo non è un errore della DTF: è precisamente la nozione di trasferimento totale descritta dalla matrice $\mathbf{H}(f)$. Diventa un errore solo quando il risultato viene etichettato come connessione anatomica o influenza diretta.

La normalizzazione per riga costruisce invece un bilancio degli ingressi del bersaglio. Se più sorgenti contribuiscono a $i$, la DTF indica la quota relativa attribuibile a ciascuna. Come per la PDC, modificare la composizione della rete può cambiare i valori normalizzati anche quando una singola connessione rimane invariata.

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

Una regola pratica è formulare verbalmente la domanda prima di scegliere la misura. “Quali bersagli ricevano un'influenza diretta da questa sorgente?” conduce naturalmente alla PDC. “Da quali sorgenti, anche attraverso percorsi intermedi, proviene l'attività che raggiunge questo bersaglio?” conduce alla DTF. Quando interessa distinguere propagazione totale e percorso diretto, le due mappe vanno lette congiuntamente invece di cercare una misura universalmente migliore.

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

L'esempio separa volutamente due ingredienti del VAR. I coefficienti diagonali di $x_1$ e $x_3$ generano oscillazioni autonome rispettivamente attorno a $f_1$ e $f_3$; i coefficienti fuori diagonale stabiliscono dove queste oscillazioni possono propagarsi. Di conseguenza un picco nato in una sorgente può comparire nello spettro di altri nodi senza che questi possiedano lo stesso oscillatore intrinseco.

Per leggere le matrici di grafici conviene procedere in tre passi: identificare prima le frequenze proprie sui pannelli diagonali; localizzare poi le PDC fuori diagonale, che ricostruiscono le frecce dirette; esaminare infine la DTF per seguire la diffusione delle oscillazioni lungo tutti i cammini. La differenza tra $1\rightarrow3$ nelle due misure rende visibile la distinzione tra struttura della rete e propagazione del segnale.

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

In questo contesto la non stazionarietà non è soltanto un disturbo da rimuovere: può essere proprio il fenomeno di interesse. Un cambiamento dei coefficienti descrive una variazione nel modo in cui il passato dei segnali contribuisce al presente. L'obiettivo passa quindi dalla stima di una singola rete media alla ricostruzione di una sequenza di reti locali.

Questa maggiore flessibilità introduce però una difficoltà di identificazione. Una variazione rapida osservata nelle mappe può derivare da un vero cambiamento della dinamica, da rumore, da un numero insufficiente di campioni o dall'algoritmo adattivo. La risoluzione temporale nominale non coincide necessariamente con quella effettiva: ogni stima contiene informazione proveniente da un intorno temporale determinato dalla finestra o dalla memoria del filtro.

### 15.1 Finestre scorrevoli

Si stima un VAR indipendente in ogni finestra. Il metodo è trasparente, ma assume stazionarietà locale. Una finestra lunga riduce la varianza e smussa le transizioni; una finestra corta migliora la localizzazione temporale ma rende la stima più instabile.

Con più prove sperimentali, una strategia comune consiste nel raccogliere, per ogni istante, i campioni appartenenti alla stessa finestra in tutte le prove. In tal modo la finestra può restare corta mentre il numero totale di osservazioni cresce. L'assunzione è che i trial siano realizzazioni dello stesso processo locale. Differenze sistematiche tra prove, jitter temporale o adattamento durante l'esperimento violano questa ipotesi e possono essere scambiati per variabilità della connettività.

Anche il passo di scorrimento merita attenzione. Finestre fortemente sovrapposte producono mappe visivamente fluide, ma le stime adiacenti non sono indipendenti. Una transizione apparentemente graduale può essere in parte una conseguenza della sovrapposizione, non della fisiologia.

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

Quando $\lambda$ è vicino a uno, il passato viene dimenticato lentamente: la stima è più regolare ma segue con ritardo i cambiamenti. Valori più piccoli aumentano la sensibilità alle osservazioni recenti e anche la variabilità. L'RLS realizza dunque in forma ricorsiva lo stesso compromesso della finestra scorrevole, sostituendo un taglio netto con una pesatura esponenziale della storia.

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

Nel TV-VAR lo stato nascosto non è l'attività neurale, ma l'insieme dei coefficienti che governano la relazione tra i segnali. Le osservazioni sono i campioni misurati e la matrice di osservazione è costruita dai loro valori passati. Il filtro aggiorna quindi, istante per istante, una distribuzione di plausibilità sui coefficienti: non produce soltanto una curva, ma anche una misura dell'incertezza associata alla curva.

La logica ricorsiva alterna due fonti di informazione. Il modello di stato esprime quanto ci si aspetta che i coefficienti persistano; il nuovo campione indica quanto la previsione basata sui coefficienti precedenti sia compatibile con i dati. Il guadagno di Kalman decide il peso relativo delle due fonti in funzione delle rispettive covarianze.

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

La predizione aumenta normalmente l'incertezza attraverso $\mathbf{Q}_{\theta}$; la correzione la riduce quando l'osservazione è informativa. Se il dato è molto rumoroso, il guadagno diminuisce e il filtro si affida maggiormente alla dinamica precedente. Se invece il modello ammette coefficienti rapidamente variabili, il filtro attribuisce più peso all'innovazione osservata.

Questa terminologia può generare confusione perché esistono due innovazioni: l'innovazione del VAR, cioè la parte del segnale non spiegata dai suoi ritardi, e l'innovazione del filtro di Kalman $\boldsymbol{\nu}_t$, cioè la differenza tra osservazione e previsione prodotta dallo stato. Sono concetti collegati ma non identici e hanno covarianze diverse.

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

La scelta del random walk non afferma che i coefficienti biologici seguano realmente un moto casuale. È una regolarizzazione dinamica: tra due istanti consecutivi ci si aspetta continuità, ma si lascia una quantità controllata di libertà per seguire i cambiamenti. Modelli di stato più strutturati sono possibili, ma richiedono ipotesi aggiuntive sulla dinamica dei coefficienti.

### 16.4 Ruolo delle covarianze

Se $\mathbf{Q}_{\theta}=q\mathbf{I}$, $q$ controlla la velocità di adattamento:

- $q$ piccolo produce coefficienti regolari ma può attenuare transizioni rapide;
- $q$ grande consente cambiamenti veloci ma può inseguire il rumore.

$\mathbf{R}_u=\operatorname{Cov}(\mathbf{u}_t)$ descrive l'incertezza dell'equazione VAR e può contenere dipendenze istantanee tra canali. Può essere fissata, inizializzata con un VAR stazionario o aggiornata ricorsivamente.

Nei dati multitrial, a ogni istante si possono usare più realizzazioni dello stesso processo, assumendo che condividano le stesse matrici $\mathbf{A}_r(t)$. La standardizzazione tra le prove deve essere descritta con precisione: media e deviazione standard vengono calcolate lungo la dimensione delle prove, separatamente per ogni tempo e canale, se questa è la procedura scelta.

Il rapporto tra le covarianze di stato e osservazione determina la scala temporale effettiva dell'adattamento. Per questo non è sufficiente riportare che è stato usato “un filtro di Kalman”: occorre descrivere inizializzazione, struttura di $\mathbf{Q}_{\theta}$, stima di $\mathbf{R}_u$ ed eventuali fattori di fading. Milde et al. (2010) mostrano inoltre che la formulazione del filtro diventa cruciale nei sistemi ad alta dimensionalità, dove una stima ricorsiva ingenua può essere numericamente o statisticamente fragile.

La taratura non dovrebbe essere scelta per rendere le mappe più nette. Va invece verificata su simulazioni con cambiamenti noti, mediante errore di previsione e attraverso analisi di sensibilità. Una connessione che appare solo per un intervallo ristretto di parametri adattivi è meno robusta di una struttura che persiste in una regione plausibile di scelte.

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

Prima di interpretare una mappa conviene separare tre livelli: variazione dei coefficienti, variazione della potenza e variazione delle misure normalizzate. Una PDC o DTF può cambiare perché cambia la connessione, ma anche perché cambiano le altre uscite o gli altri ingressi che compaiono nel denominatore. Nei dati neurofisiologici è quindi utile affiancare alle mappe direzionali i coefficienti stimati, gli spettri locali e una misura dell'incertezza.

## 17. Inferenza sulle mappe tempo-frequenza

Una mappa contiene molti test correlati. L'interpretazione di singoli pixel sopra una soglia puntuale è fragile. Quando appropriato si usano max-statistic o test basati su cluster, mantenendo la struttura spazio-temporale nei surrogati.

Prima dell'analisi conviene definire bande e finestre temporali motivate dalla domanda fisiologica. Vanno riportati ordine del VAR, metodo adattivo, parametri di regolarizzazione, numero di surrogati e famiglia di confronti controllata.

Nei dati tempo-frequenza i punti vicini sono fortemente correlati. Trattare ogni pixel come un test indipendente produce una moltiplicazione artificiale delle evidenze. I surrogati devono conservare la struttura che non appartiene all'ipotesi nulla, e la statistica massima o i cluster devono essere calcolati sull'intera famiglia di tempi, frequenze e direzioni che si intende dichiarare significativa.

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


