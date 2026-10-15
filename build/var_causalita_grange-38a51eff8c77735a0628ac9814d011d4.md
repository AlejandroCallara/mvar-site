# VAR e causalità di Granger

## 1. Dal problema biomedico alla connettività direzionale

Molti segnali biomedici sono acquisiti simultaneamente da più canali. Un EEG multicanale, per esempio, descrive l'attività elettrica osservata in più posizioni dello scalpo. La domanda che ci poniamo non riguarda soltanto la presenza di oscillazioni comuni. Spesso interessa stabilire se il passato di un processo contiene informazione utile per prevederne un altro.

Questa domanda appartiene allo studio della **connettività effettiva**, cioè delle dipendenze direzionali tra processi. Va distinta dalla connettività funzionale, che descrive associazioni statistiche senza attribuire necessariamente una direzione. Ad esempio, una correlazione elevata tra due segnali non stabilisce quale dei due anticipi l'altro e può derivare da una terza sorgente comune.

```{figure} ../assets/figures/tipi_connettivita.png
:alt: Confronto tra connettività strutturale, funzionale ed effettiva
:width: 95%
:name: fig-tipi-connettivita

Tre livelli di descrizione della connettività: struttura anatomica, dipendenza statistica e influenza direzionale. Schema tratto dalla presentazione del corso.
```

I tre livelli della figura rispondono a domande diverse. La connettività **strutturale** riguarda l'esistenza di un supporto anatomico; la connettività **funzionale** descrive il fatto che due attività varino insieme; la connettività **effettiva** cerca invece di caratterizzare un'influenza orientata. Passare da un livello all'altro non è automatico: due regioni anatomicamente connesse possono non mostrare un'interazione misurabile in una certa condizione, mentre due regioni funzionalmente associate possono non esibire una connessione anatomica diretta. La causalità di Granger fornisce una definizione operativa dell'ultimo livello, costruita sui dati osservati e non su una descrizione biofisica completa della rete (Bressler & Seth, 2011; Bastos & Schoffelen, 2016).

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

Con la convenzione adottata, $a_{21}$ descrive $1\rightarrow2$ e $a_{12}$ descrive $2\rightarrow1$. Un coefficiente non nullo indica che il passato della sorgente entra nel modello del bersaglio. Come vedremo più avanti, la causalità di Granger richiede però un confronto formale tra modelli, non la sola ispezione del coefficiente.

Per esempio, se $a_{21}\neq0$, conoscere $x_1(t-1)$ aiuta potenzialmente a prevedere $x_2(t)$. La parola “potenzialmente” è necessaria perché il coefficiente stimato può essere diverso da zero per variabilità campionaria, collinearità o specificazione inadeguata del modello. La GC valuta l'effetto congiunto di tutti i $p$ coefficienti associati alla sorgente e domanda se la loro rimozione peggiori in modo apprezzabile la previsione.

## 3. Assunzioni del modello VAR

### 3.1 Stazionarietà debole

Un processo è debolmente stazionario quando:

- la media non dipende da $t$;
- la varianza è finita e costante;
- la covarianza tra $\mathbf{x}(t)$ e $\mathbf{x}(t-h)$ dipende dal ritardo $h$, non dall'istante $t$.

Trend, cambiamenti di varianza e transizioni tra stati fisiologici violano questa ipotesi. Prima della stima possono servire detrending, rimozione della media, trasformazioni o una formulazione tempo-variante.

La stazionarietà deve essere valutata rispetto all’intervallo sul quale viene stimato il modello. Una registrazione EEG può cambiare nel corso dell’esperimento, ma risultare approssimativamente stazionaria all’interno di una breve fase sperimentale. Accorciare la finestra consente di seguire meglio questi cambiamenti, ma riduce il numero di campioni disponibili e aumenta l’incertezza della stima VAR. La durata della finestra deve quindi bilanciare omogeneità della dinamica, precisione temporale e affidabilità statistica.

È inoltre utile distinguere la stazionarietà dei dati dalla stabilità del modello stimato. La stazionarietà riguarda il segmento osservato e richiede che le sue proprietà statistiche, come media, varianza e correlazioni temporali, rimangano approssimativamente costanti. La stabilità riguarda invece la dinamica del modello stimato: gli effetti di una perturbazione devono attenuarsi nel tempo e non produrre una risposta divergente. Un VAR può risultare stabile anche quando è stato stimato su un segmento contenente cambiamenti di stato, trend o variazioni di potenza. La verifica delle radici del modello deve quindi essere accompagnata dall’esame dei segnali, delle statistiche calcolate in intervalli successivi e dei residui.

### 3.2 Stabilità

La stabilità è una proprietà del modello. Un VAR stabile dissipa l'effetto delle condizioni iniziali e delle perturbazioni. Data la matrice compagna associata al VAR:

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

Nel modello VAR, il vettore delle innovazioni $\mathbf{u}(t)$ rappresenta la parte dei segnali al tempo $t$ che non può essere prevista utilizzando il passato delle variabili incluse nel modello. Dopo avere calcolato la previsione

$$
\widehat{\mathbf{x}}(t)
=
\mathbf{c}
+
\sum_{r=1}^{p}\mathbf{A}_r\mathbf{x}(t-r),
$$

l'innovazione è definita come l'errore di previsione a un passo:

$$
\mathbf{u}(t)
=
\mathbf{x}(t)-\widehat{\mathbf{x}}(t).
$$

Il termine *innovazione* sottolinea che questa componente contiene informazione nuova rispetto al passato considerato dal modello. Non coincide necessariamente con il rumore dello strumento: può comprendere attività fisiologica non prevedibile, ingressi esterni non misurati, variabili omesse, dinamiche non lineari e rumore di misura.

Se il modello descrive adeguatamente la dinamica temporale, le innovazioni non devono conservare una struttura prevedibile nel tempo. In particolare, per $h\neq0$ dovrebbe valere approssimativamente

$$
\operatorname{Cov}
\left[
\mathbf{u}(t),\mathbf{u}(t-h)
\right]
=
\mathbf{0}.
$$

Innovazioni serialmente correlate indicano che nei residui è rimasta informazione temporale. Ciò può suggerire, per esempio, un ordine VAR insufficiente, la presenza di variabili mancanti o una dinamica non adeguatamente rappresentata dal modello.

Le innovazioni dei diversi canali possono invece essere correlate nello stesso istante. La loro covarianza è descritta dalla matrice

$$
\boldsymbol{\Sigma}_u
=
\operatorname{Cov}
\left[
\mathbf{u}(t)
\right].
$$

Un elemento fuori diagonale $\Sigma_{u,ij}\neq0$ indica che le equazioni dei canali $i$ e $j$ ricevono componenti imprevedibili contemporaneamente correlate. Questa dipendenza può derivare da un ingresso comune non osservato, da interazioni più rapide della risoluzione temporale o dal mixing dei segnali.

È importante distinguere questa covarianza istantanea dalle influenze ritardate descritte dalle matrici $\mathbf{A}_r$. Poiché le due innovazioni si manifestano nello stesso campione, il VAR non dispone dell'ordine temporale necessario per stabilire se l'effetto vada da $i$ a $j$ oppure da $j$ a $i$. Una covarianza residua fuori diagonale segnala quindi una dipendenza istantanea, ma non identifica da sola una direzione causale (Faes & Nollo, 2011).

Infine, si assume che le innovazioni abbiano media nulla e covarianza finita. L'ipotesi di gaussianità non è indispensabile per interpretare la stima ai minimi quadrati, ma rende più semplice formulare la verosimiglianza e applicare alcuni test statistici parametrici.

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

Con dati multitrial (come ad esempio gli ERP), i campioni validi delle diverse prove possono contribuire alla stessa regressione, purché le prove siano considerate realizzazioni comparabili dello stesso processo. Non si devono però creare ritardi artificiali concatenando la fine di una prova con l'inizio della successiva: i regressori vanno costruiti separatamente in ogni trial e solo dopo impilati. Questa distinzione diventa decisiva nelle finestre brevi, dove l'uso di più prove fornisce le osservazioni necessarie senza allungare la finestra temporale.

## 5. Scelta dell'ordine del modello

L'ordine $p$ indica quanti campioni passati vengono utilizzati dal modello VAR per prevedere il campione corrente. Determina quindi la memoria temporale del modello:

$$
\tau_{\mathrm{VAR}}=\frac{p}{F_s},
$$

dove $F_s$ è la frequenza di campionamento. Per esempio, con $F_s=250$ Hz e $p=10$, il modello utilizza una storia temporale di

$$
\tau_{\mathrm{VAR}}=\frac{10}{250}=0.04\ \mathrm{s}=40\ \mathrm{ms}.
$$

La scelta di $p$ deve bilanciare capacità descrittiva e complessità:

- se $p$ è troppo basso, il modello non rappresenta tutta la dipendenza temporale presente nei dati e lascia una struttura prevedibile nei residui;
- se $p$ è troppo alto, il numero di parametri aumenta, le stime diventano più variabili e il modello può adattarsi anche alle fluttuazioni casuali del campione.

```{figure} ../assets/figures/ordine_modello.png
:alt: Compromesso nella scelta dell'ordine del modello autoregressivo
:width: 80%
:name: fig-ordine-modello

Un ordine troppo basso lascia struttura temporale nei residui; un ordine troppo alto aumenta la complessità, la varianza delle stime e il rischio di sovradattamento. Schema tratto dalla presentazione del corso.
```

### 5.1 Adattamento e complessità

Per ogni ordine candidato $p$ si stima un nuovo modello VAR e si calcola la matrice di covarianza delle innovazioni

$$
\widehat{\boldsymbol{\Sigma}}_u(p).
$$

Il termine

$$
\log\det\widehat{\boldsymbol{\Sigma}}_u(p)
$$

misura l'errore di previsione multivariato. Aumentando l'ordine, questo errore tende a diminuire perché il modello dispone di un numero maggiore di coefficienti.

La riduzione dell'errore, tuttavia, non implica necessariamente che il modello descriva meglio il processo reale. Un modello più complesso può adattarsi anche al rumore. I criteri di informazione aggiungono quindi all'errore di previsione una penalizzazione che cresce con il numero di parametri:

$$
IC(p)
=
\text{errore di previsione}
+
\text{penalizzazione della complessità}.
$$

Per un VAR con $M$ variabili e ordine $p$, il numero dei coefficienti dinamici è

$$
M^2p.
$$

Ogni nuovo ritardo aggiunge infatti una matrice $\mathbf{A}_r$ contenente $M^2$ coefficienti.

### 5.2 Criteri di informazione più comuni

Indichiamo con:

- $N$ il numero di osservazioni effettivamente utilizzate per la stima;
- $M$ il numero di variabili;
- $p$ l'ordine del modello;
- $\widehat{\boldsymbol{\Sigma}}_u(p)$ la covarianza delle innovazioni ottenuta con l'ordine $p$.

| Criterio | Formula | Interpretazione |
|---|---|---|
| Akaike Information Criterion (AIC) | $\displaystyle \log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{2}{N}pM^2$ | Applica una penalizzazione moderata e può selezionare ordini relativamente elevati. |
| Schwarz-Bayes Criterion (SBC o BIC) | $\displaystyle \log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{\log N}{N}pM^2$ | Applica una penalizzazione più forte e tende a scegliere modelli più parsimoniosi. |
| Hannan-Quinn Criterion (HQ) | $\displaystyle \log\det\widehat{\boldsymbol{\Sigma}}_u(p)+\frac{2\log\log N}{N}pM^2$ | Presenta generalmente un comportamento intermedio tra AIC e BIC. |
| Akaike's Final Prediction Error (FPE) | $\displaystyle \left(\frac{N+Mp+1}{N-Mp-1}\right)^M\det\widehat{\boldsymbol{\Sigma}}_u(p)$ | Stima l'errore di predizione finale ed è strettamente collegato all'AIC. |

Per ogni criterio si sceglie l'ordine corrispondente al valore minimo:

$$
\widehat{p}
=
\underset{p}{\operatorname{arg\,min}}\;IC(p).
$$

Il minimo identifica il punto nel quale l'ulteriore riduzione dell'errore di previsione non compensa più l'aumento della complessità.

### 5.3 Interpretazione pratica

I diversi criteri possono indicare ordini differenti perché applicano penalizzazioni diverse. L'AIC privilegia maggiormente la capacità predittiva e tende spesso a scegliere ordini più elevati. Il BIC penalizza maggiormente il numero di parametri e favorisce modelli più semplici. HQ fornisce generalmente una scelta intermedia, mentre FPE è strettamente collegato all'AIC.

Il risultato di un criterio di informazione non deve essere accettato automaticamente. Se il valore continua a diminuire fino al massimo ordine esaminato, il minimo non è stato realmente individuato. In questo caso occorre ampliare l'intervallo di ricerca oppure verificare la presenza di problemi di preprocessing, non stazionarietà o inadeguatezza del modello.

Quando la curva non presenta un minimo netto, il criterio del gomito può suggerire l'ordine oltre il quale il miglioramento diventa modesto. Si tratta però di una scelta euristica che deve essere confermata dalla validazione del modello.

Anche la scala temporale del fenomeno fornisce un'indicazione utile. Se il massimo ritardo fisiologico che si desidera rappresentare è $\tau_{\max}$, l'intervallo di ordini esaminato dovrebbe consentire almeno

$$
\frac{p}{F_s}\geq\tau_{\max}.
$$

Questa relazione aiuta a definire l'intervallo di ricerca, ma non determina da sola l'ordine corretto. Un ordine abbastanza elevato da coprire il ritardo atteso potrebbe comunque produrre stime instabili o richiedere più dati di quelli disponibili.

### 5.4 Procedura consigliata

Una procedura pratica consiste nel:

1. definire un intervallo di ordini compatibile con la frequenza di campionamento e con i ritardi fisiologici di interesse;
2. stimare un modello distinto per ogni ordine candidato;
3. calcolare AIC, BIC, HQ o FPE;
4. individuare i minimi e confrontare gli ordini suggeriti;
5. validare i modelli candidati controllando stabilità e bianchezza dei residui.

Un ordine selezionato da un criterio di informazione rappresenta quindi un candidato, non una garanzia di adeguatezza. Se un modello di ordine più basso produce già residui bianchi e risultati stabili, può essere preferibile a un modello più complesso indicato dall'AIC. Al contrario, un ordine scelto dal BIC non è adeguato se lascia una chiara dipendenza temporale nei residui.

La scelta dell'ordine influenza anche l'analisi spettrale. Un ordine troppo basso può non rappresentare oscillazioni strette o ritardi distribuiti. Un ordine troppo alto può introdurre picchi instabili e relazioni direzionali spurie. La scelta finale deve quindi considerare insieme criteri di informazione, validazione statistica e interpretabilità fisiologica.


## 6. Validazione del modello

La scelta dell'ordine mediante AIC, BIC, HQ o FPE individua un modello candidato, ma non dimostra che il modello descriva adeguatamente i dati. Dopo la stima occorre verificare se il VAR ha catturato la struttura temporale rilevante e se la dinamica ottenuta è matematicamente ammissibile.

La validazione risponde quindi a tre domande:

1. nei residui è rimasta informazione temporale prevedibile?
2. il modello riproduce la struttura di correlazione osservata nei dati?
3. la dinamica stimata è stabile?

Queste domande sono complementari. Un modello stabile può lasciare autocorrelazione nei residui; un modello con residui compatibili con rumore bianco può riprodurre in modo incompleto alcune caratteristiche dei dati. Nessun controllo, considerato isolatamente, certifica la validità del modello.

### 6.1 Bianchezza dei residui: test di Portmanteau

I residui del modello stimato sono

$$
\widehat{\mathbf{u}}(t)
=
\mathbf{x}(t)-\widehat{\mathbf{x}}(t).
$$

Essi rappresentano la parte del campione corrente che il VAR non riesce a prevedere utilizzando il passato delle variabili incluse. Se il modello ha descritto adeguatamente la dipendenza temporale, i residui non dovrebbero contenere ulteriore informazione prevedibile.

I test di Portmanteau verificano congiuntamente la presenza di autocorrelazione residua fino a un ritardo massimo $h$. Le ipotesi sono

$$
H_0:
\mathbf{C}_1=\mathbf{C}_2=\cdots=\mathbf{C}_h=\mathbf{0},
$$

$$
H_1:
\mathbf{C}_{\ell}\neq\mathbf{0}
\quad
\text{per almeno un } \ell\leq h,
$$

dove $\mathbf{C}_{\ell}$ è la matrice di autocovarianza dei residui al ritardo $\ell$:

$$
\mathbf{C}_{\ell}
=
\frac{1}{T}
\sum_{t=\ell+1}^{T}
\widehat{\mathbf{u}}(t)
\widehat{\mathbf{u}}(t-\ell)^{\mathsf T}.
$$

Per riassumere la dipendenza residua al ritardo $\ell$, definiamo

$$
S_{\ell}
=
\operatorname{tr}
\left(
\mathbf{C}_{\ell}^{\mathsf T}
\mathbf{C}_0^{-1}
\mathbf{C}_{\ell}
\mathbf{C}_0^{-1}
\right).
$$

I test più comuni combinano i valori $S_{\ell}$ in una singola statistica:

| Test | Statistica | Caratteristiche |
|---|---|---|
| Box-Pierce | $\displaystyle Q_h^{\mathrm{BP}}=T\sum_{\ell=1}^{h}S_{\ell}$ | È la formulazione originale. È semplice, ma può avere proprietà poco soddisfacenti con campioni piccoli. |
| Ljung-Box | $\displaystyle Q_h^{\mathrm{LB}}=T(T+2)\sum_{\ell=1}^{h}\frac{S_{\ell}}{T-\ell}$ | Corregge la statistica di Box-Pierce per migliorare il comportamento nei piccoli campioni. |
| Li-McLeod | $\displaystyle Q_h^{\mathrm{LM}}=T\sum_{\ell=1}^{h}S_{\ell}+\frac{M^2h(h+1)}{2T}$ | Introduce una correzione specifica per il caso multivariato e presenta generalmente un comportamento migliore nei piccoli campioni. |

Sotto $H_0$ e per campioni sufficientemente grandi, la statistica viene confrontata approssimativamente con una distribuzione chi-quadrato:

$$
Q_h
\overset{a}{\sim}
\chi^2_{M^2(h-p)},
\qquad h>p,
$$

dove $M$ è il numero di variabili e $p$ è l'ordine del VAR.

L'interpretazione del risultato è:

- se il $p$-value è minore o uguale al livello di significatività $\alpha$, si rifiuta $H_0$: nei residui rimane dipendenza temporale e il modello deve essere rivisto;
- se il $p$-value è maggiore di $\alpha$, non si rifiuta $H_0$: i residui sono compatibili con l'ipotesi di rumore bianco fino al ritardo $h$.

Un $p$-value elevato non dimostra che i residui siano bianchi. Indica soltanto che il test non ha rilevato autocorrelazioni residue fino al ritardo considerato.

La scelta di $h$ influenza il risultato. Un valore troppo basso può non rilevare dipendenze a lungo termine; un valore molto elevato riduce l'affidabilità dell'approssimazione e può diminuire la potenza del test. È quindi opportuno dichiarare il valore di $h$ e verificare che la conclusione rimanga simile per alcuni valori ragionevoli.

Con dati multitrial, le autocovarianze residue devono essere calcolate rispettando i confini delle prove. Non bisogna considerare consecutivi il residuo finale di un trial e quello iniziale del trial successivo.

### 6.2 Consistenza del modello

La bianchezza dei residui verifica se nel segnale è rimasta dipendenza temporale prevedibile. La consistenza considera invece quanto bene il VAR riproduca la struttura di auto e cross-correlazione osservata nei dati.

La procedura consiste nel:

1. simulare nuovi dati utilizzando i coefficienti e la covarianza delle innovazioni del VAR stimato;
2. calcolare le auto e cross-correlazioni dei dati reali e di quelli simulati;
3. confrontare le due strutture mediante un indice di somiglianza.

Indicando con $\mathbf{R}_{\mathrm{real}}$ e $\mathbf{R}_{\mathrm{sim}}$ i vettori contenenti le correlazioni dei dati reali e simulati, la percent consistency può essere definita come

$$
PC
=
100
\left(
1-
\frac{
\left\|
\mathbf{R}_{\mathrm{sim}}-\mathbf{R}_{\mathrm{real}}
\right\|_2
}{
\left\|
\mathbf{R}_{\mathrm{real}}
\right\|_2
}
\right).
$$

Un valore elevato indica che il modello genera dati con una struttura correlativa simile a quella osservata. Un valore basso suggerisce che il VAR non riproduce adeguatamente alcune dipendenze presenti nei dati.

La consistenza dipende dai ritardi utilizzati nel confronto e dalla variabilità della simulazione. Non esiste quindi una soglia universale valida per ogni applicazione. Inoltre, una buona consistenza globale non garantisce che ogni singola relazione direzionale sia corretta.

### 6.3 Stabilità del modello

La stabilità riguarda il comportamento dinamico del VAR. Un modello stabile dissipa nel tempo gli effetti delle condizioni iniziali e delle perturbazioni. Si verifica calcolando gli autovalori della matrice compagna stimata:

$$
\rho\!\left(\widehat{\mathbf{A}}_c\right)
=
\max_k
\left|
\lambda_k
\right|
<1.
$$

Se un autovalore si trova sul cerchio unitario o al suo esterno, il modello non è stabile e può generare risposte persistenti o divergenti. In questo caso le misure di connettività derivate dal VAR non dovrebbero essere interpretate.

La stabilità non garantisce però che il modello sia adeguato. Anche un VAR stabile può avere un ordine insufficiente, lasciare dipendenza nei residui o descrivere male la struttura di correlazione dei dati.

### 6.4 Valutazione complessiva

Un modello può essere considerato un candidato adeguato quando:

- è stabile;
- i test di Portmanteau non rilevano una dipendenza temporale residua significativa;
- riproduce in modo soddisfacente la struttura correlativa dei dati;
- i risultati rimangono ragionevolmente simili per scelte vicine dell'ordine del modello.

Se uno di questi controlli fallisce, occorre riesaminare l'ordine, l'insieme delle variabili, la stazionarietà del segmento e le scelte di preprocessing prima di interpretare le misure di causalità.

La [Demo 45](../demo/demo_45_validazione.md) mostra perché la sola stabilità non è sufficiente: un VAR con ordine sottostimato può essere stabile ma fallire i test sui residui.

## 7. Causalità di Granger nel dominio del tempo

La causalità di Granger traduce una domanda sulla direzione dell'informazione in un confronto tra previsioni:

> Il passato di una variabile migliora la previsione di un'altra variabile, una volta considerata l'informazione già contenuta nel resto del sistema?

Supponiamo di voler verificare la direzione $x_j\rightarrow x_i$. Si confrontano due modelli del bersaglio $x_i$:

- un modello che utilizza anche il passato della sorgente $x_j$;
- un modello dal quale il passato di $x_j$ è stato escluso.

Se l'esclusione peggiora la previsione, il passato di $x_j$ contiene informazione utile e non ridondante per prevedere $x_i$. In questo caso si dice che $x_j$ causa $x_i$ in senso di Granger.

Questa definizione riguarda la predicibilità dei segnali osservati. Da sola non dimostra un meccanismo fisiologico diretto e non esclude l'effetto di variabili non misurate.

### 7.1 Modello completo e modello ridotto

Consideriamo inizialmente due processi $x_1$ e $x_2$ e vogliamo verificare la direzione

$$
x_1\rightarrow x_2.
$$

Il **modello completo** utilizza sia il passato di $x_2$ sia quello di $x_1$:

$$
x_2(t)
=
\sum_{r=1}^{p}a_{22,r}x_2(t-r)
+
\sum_{r=1}^{p}a_{21,r}x_1(t-r)
+
u_{2,F}(t).
$$

La varianza dell'errore di previsione del modello completo è

$$
\sigma_{2,F}^{2}
=
\operatorname{var}\!\left[u_{2,F}(t)\right].
$$

Il **modello ridotto** esclude tutti i ritardi di $x_1$:

$$
x_2(t)
=
\sum_{r=1}^{p}\widetilde{a}_{22,r}x_2(t-r)
+
u_{2,R}(t),
$$

con varianza residua

$$
\sigma_{2,R}^{2}
=
\operatorname{var}\!\left[u_{2,R}(t)\right].
$$

I coefficienti del modello ridotto devono essere ristimati. Non è sufficiente stimare il modello completo e porre successivamente a zero i coefficienti $a_{21,r}$, perché l'esclusione della sorgente modifica anche la stima dei coefficienti associati agli altri regressori.

L'ipotesi nulla è

$$
H_0:
a_{21,1}=a_{21,2}=\cdots=a_{21,p}=0.
$$

Sotto $H_0$, il passato di $x_1$ non fornisce informazione aggiuntiva per prevedere $x_2$. L'ipotesi alternativa richiede che almeno uno di questi coefficienti sia diverso da zero.

La direzione opposta deve essere verificata separatamente. Per studiare $x_2\rightarrow x_1$ si escludono tutti i ritardi di $x_2$ dall'equazione di $x_1$.

### 7.2 Misura della causalità di Granger

La causalità di Granger nel dominio del tempo quantifica il peggioramento della previsione mediante il rapporto tra le varianze residue:

$$
F_{1\rightarrow2}
=
\log
\frac{\sigma_{2,R}^{2}}
{\sigma_{2,F}^{2}}.
$$

L'interpretazione è immediata:

- se $\sigma_{2,R}^{2}\approx\sigma_{2,F}^{2}$, l'esclusione di $x_1$ non modifica la previsione e $F_{1\rightarrow2}\approx0$;
- se $\sigma_{2,R}^{2}>\sigma_{2,F}^{2}$, il modello completo predice meglio $x_2$ e $F_{1\rightarrow2}>0$;
- quanto maggiore è l'aumento della varianza residua, tanto maggiore è il valore della GC.

Per esempio, se

$$
\sigma_{2,F}^{2}=0.8,
\qquad
\sigma_{2,R}^{2}=1,
$$

allora

$$
F_{1\rightarrow2}
=
\log\left(\frac{1}{0.8}\right)
\approx0.223.
$$

In questo caso l'esclusione del passato di $x_1$ aumenta del 25% la varianza dell'errore di previsione:

$$
\exp(F_{1\rightarrow2})
=
\frac{\sigma_{2,R}^{2}}{\sigma_{2,F}^{2}}
=1.25.
$$

Con modelli annidati, stimati sugli stessi campioni e con lo stesso metodo, ci si attende

$$
\sigma_{2,R}^{2}\geq\sigma_{2,F}^{2}
$$

e quindi una GC non negativa. Piccoli valori negativi possono comparire per errori numerici, regolarizzazione o procedure di stima non perfettamente confrontabili.

La GC non coincide con il valore di un singolo coefficiente. Il test riguarda congiuntamente tutti i coefficienti associati ai $p$ ritardi della sorgente. Inoltre, il valore della misura dipende dalla dinamica del bersaglio, dalla covarianza delle innovazioni e dall'informazione condivisa con le altre variabili.

### 7.3 Causalità di Granger condizionata

Manteniamo la notazione utilizzata nelle slide. Indichiamo con

$$
V=\{1,2,\ldots,M\}
$$

l'insieme di tutte le variabili osservate. Vogliamo verificare se $x_j$ Granger-causa $x_i$, condizionatamente alle variabili appartenenti a $V\setminus\{i,j\}$.

Il **modello completo** include il passato di tutte le variabili:

$$
\mathbf{x}_t^V
=
\sum_{k=1}^{p}
\mathbf{A}_k
\mathbf{x}_{t-k}^V
+
\mathbf{u}_t.
$$

La varianza dell'errore nella previsione di $x_i$ è l'elemento diagonale corrispondente della matrice di covarianza delle innovazioni:

$$
\operatorname{var}(u_{i,t})
=
\Sigma_{ii}.
$$

Il **modello ridotto** esclude $x_j$ dall'insieme dei predittori:

$$
\mathbf{x}_t^{V\setminus j}
=
\sum_{k=1}^{p}
\overline{\mathbf{A}}_k
\mathbf{x}_{t-k}^{V\setminus j}
+
\overline{\mathbf{u}}_t.
$$

La barra distingue i parametri e le innovazioni del modello ridotto da quelli del modello completo. Tutti i coefficienti del modello ridotto vengono ristimati. La varianza dell'errore nella previsione di $x_i$ diventa

$$
\operatorname{var}(\overline{u}_{i,t})
=
\overline{\Sigma}_{ii}.
$$

Poiché il modello ridotto dispone di meno informazioni, ci si attende

$$
\overline{\Sigma}_{ii}\geq\Sigma_{ii}.
$$

La causalità di Granger condizionata è quindi definita come

$$
F_{j\rightarrow i\mid V\setminus\{i,j\}}
=
\ln
\frac{\overline{\Sigma}_{ii}}
{\Sigma_{ii}}.
$$

Un valore

$$
F_{j\rightarrow i\mid V\setminus\{i,j\}}>0
$$

indica che il passato di $x_j$ migliora la previsione di $x_i$, oltre alle informazioni già contenute nel passato delle altre variabili.

In modo equivalente, l'assenza di causalità di Granger da $x_j$ a $x_i$ corrisponde a

$$
x_j\nrightarrow x_i
\quad\Longleftrightarrow\quad
A_{k,ij}=0
\qquad
\forall k=1,\ldots,p.
$$

Pertanto, $x_j$ Granger-causa $x_i$ quando almeno uno dei coefficienti $A_{k,ij}$, con $k=1,\ldots,p$, risulta significativamente diverso da zero.

Il condizionamento permette di distinguere alcune relazioni dirette da effetti indiretti. Consideriamo, per esempio, il percorso

$$
x_1\rightarrow x_3\rightarrow x_2.
$$

Un modello bivariato contenente soltanto $x_1$ e $x_2$ può rilevare $x_1\rightarrow x_2$. Inserendo $x_3$, la relazione può ridursi o scomparire perché l'informazione trasmessa da $x_1$ raggiunge $x_2$ attraverso $x_3$.

Questo vantaggio dipende però dalle variabili osservate. Se una sorgente comune non viene misurata, il condizionamento non può eliminarne l'effetto. Anche l'inclusione di variabili molto rumorose o scarsamente motivate può rendere la stima meno affidabile. L'insieme dei segnali deve quindi essere scelto sulla base del disegno sperimentale e delle ipotesi fisiologiche.

### 7.4 Dipendenza istantanea

Le matrici dei coefficienti VAR descrivono dipendenze ritardate. Le innovazioni di due equazioni possono però essere correlate nello stesso campione:

$$
\Sigma_{u,ij}
=
\operatorname{Cov}
\left[
u_i(t),u_j(t)
\right]
\neq0.
$$

Questa covarianza indica che i due segnali ricevono componenti nuove contemporaneamente correlate. Può derivare da:

- un ingresso comune non osservato;
- un'interazione più rapida della risoluzione temporale disponibile;
- mixing istantaneo o conduzione di volume;
- rumore di misura correlato.

Nel caso bivariato, una misura della dipendenza istantanea è

$$
F_{i\cdot j}
=
\log
\frac{
\Sigma_{u,ii}\Sigma_{u,jj}
}{
\det
\begin{bmatrix}
\Sigma_{u,ii} & \Sigma_{u,ij}\\
\Sigma_{u,ji} & \Sigma_{u,jj}
\end{bmatrix}
}.
$$

Questa misura è simmetrica:

$$
F_{i\cdot j}=F_{j\cdot i}.
$$

La dipendenza istantanea segnala una relazione tra le innovazioni, ma non permette di assegnare una direzione. Poiché le due variazioni avvengono nello stesso campione, manca una precedenza temporale che consenta di distinguere $i\rightarrow j$ da $j\rightarrow i$.

### 7.5 Inferenza parametrica

La misura

$$
F_{j\rightarrow i\mid V\setminus\{i,j\}}
=
\ln
\frac{\overline{\Sigma}_{ii}}
{\Sigma_{ii}}
$$

è teoricamente non negativa. Occorre però stabilire se la riduzione della varianza residua sia abbastanza grande da indicare un contributo predittivo reale oppure se possa essere spiegata dalla variabilità campionaria.

L'inferenza parametrica viene formulata mediante un **likelihood-ratio test**. L'ipotesi nulla è l'assenza di causalità di Granger:

$$
H_0:
F_{j\rightarrow i\mid V\setminus\{i,j\}}=0.
$$

Questa ipotesi equivale a imporre che tutti i coefficienti che collegano i ritardi di $x_j$ a $x_i$ siano nulli:

$$
H_0:
A_{1,ij}=A_{2,ij}=\cdots=A_{p,ij}=0.
$$

Il modello ridotto impone quindi $p$ vincoli al modello completo. Sotto $H_0$ e per un numero sufficientemente grande di osservazioni, la statistica del likelihood ratio segue approssimativamente una distribuzione chi-quadrato:

$$
N
F_{j\rightarrow i\mid V\setminus\{i,j\}}
\overset{a}{\sim}
\chi_p^2,
$$

dove $N$ è il numero di osservazioni effettivamente utilizzate per la stima. Per una singola serie di lunghezza $T$,

$$
N=T-p.
$$

Il numero di gradi di libertà è $p$ perché l'ipotesi nulla impone un vincolo per ciascun ritardo della sorgente.

Il $p$-value si ottiene confrontando il valore osservato

$$
N
F_{j\rightarrow i\mid V\setminus\{i,j\}}^{\mathrm{obs}}
$$

con la coda superiore della distribuzione $\chi_p^2$. Un $p$-value piccolo fornisce evidenza contro $H_0$: il modello completo offre un miglioramento predittivo maggiore di quello atteso per effetto della variabilità campionaria.

Questa inferenza è asintotica. Con pochi campioni, residui fortemente non gaussiani o un modello non adeguatamente validato, l'approssimazione chi-quadrato può risultare poco affidabile. In questi casi è utile ricorrere a una distribuzione nulla empirica, come quella ottenuta mediante bootstrap.

La significatività statistica deve essere distinta dall'entità dell'effetto. Con molti campioni anche un valore molto piccolo di GC può risultare significativo; con pochi campioni un effetto rilevante può non superare la soglia statistica. È quindi opportuno riportare sia il valore della GC sia il risultato del test.

Quando si analizzano molte direzioni, finestre temporali o frequenze, occorre inoltre correggere per confronti multipli.

### 7.6 Bootstrap sotto l'ipotesi nulla

Il bootstrap permette di costruire empiricamente la distribuzione della GC attesa in assenza della connessione. Per verificare $x_j\rightarrow x_i$ si procede nel modo seguente:

1. si stima il modello ridotto, nel quale i ritardi di $x_j$ sono esclusi dall'equazione di $x_i$;
2. si ricampionano i vettori delle innovazioni del modello ridotto;
3. si simulano nuovi dati sotto l'ipotesi nulla;
4. per ogni dataset simulato si ristimano sia il modello completo sia quello ridotto;
5. si calcola la GC ottenendo una distribuzione empirica sotto $H_0$;
6. si confronta la GC osservata con questa distribuzione.

I residui devono essere ricampionati come vettori multivariati, non separatamente per ogni canale. In questo modo si conserva la covarianza istantanea tra le innovazioni.

Con $B$ repliche, il $p$-value empirico può essere calcolato come

$$
p_{\mathrm{boot}}
=
\frac{
1+\#\left\{
F_b\geq F_{\mathrm{obs}}
\right\}
}{
B+1
}.
$$

L'aggiunta di 1 al numeratore e al denominatore evita di ottenere un $p$-value esattamente uguale a zero quando nessuna replica supera il valore osservato.

La [Demo 50](../demo/demo_50_bootstrap.md) applica questa procedura alle due possibili direzioni di un sistema bivariato.

### 7.7 Interpretazione

Una relazione $x_j\rightarrow x_i$ indica che il passato di $x_j$ contiene informazione utile per prevedere $x_i$, condizionatamente alle variabili incluse nel modello. Non dimostra da sola una connessione anatomica o un meccanismo fisiologico diretto.

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

Molti segnali fisiologici presentano attività oscillatoria. Nell'EEG, per esempio, alcune interazioni possono manifestarsi soprattutto in specifiche bande di frequenza e risultare poco evidenti quando vengono riassunte da un'unica misura nel dominio del tempo.

La rappresentazione spettrale permette di rispondere a domande come:

- a quali frequenze due segnali mostrano una relazione?
- la relazione è direzionale?
- il contributo osservato rappresenta una connessione diretta o la propagazione attraverso nodi intermedi?
- il flusso viene interpretato rispetto alle uscite della sorgente o rispetto agli ingressi del bersaglio?

L'analisi in frequenza non richiede di stimare un nuovo modello. Si stima e si valida il VAR nel dominio del tempo, poi si trasformano i coefficienti stimati per descrivere la stessa dinamica frequenza per frequenza. Le misure spettrali ereditano quindi tutte le ipotesi e gli eventuali limiti del VAR di partenza. Un modello instabile, residui non bianchi o un ordine inadeguato compromettono anche PDC, DTF e le altre quantità derivate.

### 9.1 Dal modello temporale a $\mathbf{A}(f)$

Consideriamo il modello VAR($p$), dopo avere rimosso la media:

$$
\mathbf{x}_t
=
\sum_{k=1}^{p}
\mathbf{A}_k\mathbf{x}_{t-k}
+
\mathbf{u}_t.
$$

Riscrivendolo in modo da isolare le innovazioni:

$$
\mathbf{u}_t
=
\mathbf{x}_t
-
\sum_{k=1}^{p}
\mathbf{A}_k\mathbf{x}_{t-k}.
$$

Applicando la trasformata di Fourier si ottiene

$$
\mathbf{U}(f)
=
\mathbf{A}(f)\mathbf{X}(f),
$$

dove

$$
\mathbf{A}(f)
=
\mathbf{I}
-
\sum_{k=1}^{p}
\mathbf{A}_k
e^{-i2\pi fk}.
$$

Questa è la notazione adottata nelle slide, nella quale $f$ è espresso come frequenza normalizzata. Se si usa la frequenza fisica in hertz e il periodo di campionamento è $T_s$, il termine esponenziale diventa $e^{-i2\pi f kT_s}$.

La matrice $\mathbf{A}(f)$ contiene la rappresentazione in frequenza dei coefficienti autoregressivi. L'elemento $A_{ij}(f)$ descrive il legame diretto dalla sorgente $j$ al bersaglio $i$ alle diverse frequenze, condizionatamente alle altre variabili incluse nel modello.

### 9.2 Matrice di trasferimento

Premoltiplicando per $\mathbf{A}(f)^{-1}$ si ottiene

$$
\mathbf{X}(f)
=
\mathbf{A}(f)^{-1}\mathbf{U}(f)
=
\mathbf{H}(f)\mathbf{U}(f),
$$

dove

$$
\mathbf{H}(f)
=
\mathbf{A}(f)^{-1}
$$

è la matrice di trasferimento.

La matrice $\mathbf{H}(f)$ descrive come le innovazioni si propagano attraverso il sistema. L'elemento $H_{ij}(f)$ quantifica il contributo dell'innovazione generata nel nodo $j$ al segnale osservato nel nodo $i$. Poiché l'inversione di $\mathbf{A}(f)$ incorpora tutti i cammini presenti nella rete, questo contributo può comprendere sia percorsi diretti sia percorsi indiretti.

Le due matrici rispondono quindi a domande differenti:

- $\mathbf{A}(f)$ descrive i legami diretti imposti dai coefficienti del VAR;
- $\mathbf{H}(f)$ descrive la propagazione complessiva delle innovazioni attraverso la rete.

Questa distinzione costituisce la base della differenza tra PDC e DTF.

### 9.3 Matrice di densità spettrale

Indichiamo con

$$
\boldsymbol{\Sigma}
=
\operatorname{Cov}(\mathbf{u}_t)
$$

la matrice di covarianza delle innovazioni. La matrice di densità spettrale del processo è

$$
\mathbf{S}(f)
=
\operatorname{E}
\left[
\mathbf{X}(f)\mathbf{X}(f)^{\mathrm H}
\right]
=
\mathbf{H}(f)
\boldsymbol{\Sigma}
\mathbf{H}(f)^{\mathrm H},
$$

dove $\mathrm{H}$ indica il trasposto coniugato.

Gli elementi di $\mathbf{S}(f)$ hanno due ruoli:

- gli elementi diagonali $S_{ii}(f)$ sono gli auto-spettri e descrivono la potenza della variabile $i$ alla frequenza $f$;
- gli elementi fuori diagonale $S_{ij}(f)$ sono i cross-spettri e descrivono ampiezza e fase della relazione tra le variabili $i$ e $j$.

A partire da $\mathbf{S}(f)$, $\mathbf{A}(f)$, $\mathbf{H}(f)$ e $\boldsymbol{\Sigma}$ si possono costruire misure con significati differenti:

- spettro e cross-spettro descrivono potenza e relazioni spettrali;
- coerenza e coerenza parziale descrivono associazioni non direzionali;
- PDC descrive l'influenza diretta e condizionata associata a $\mathbf{A}(f)$;
- DTF descrive la propagazione complessiva associata a $\mathbf{H}(f)$.

La frequenza alla quale si valuta una misura non crea nuova informazione rispetto al VAR. Essa mostra come la dinamica descritta dai coefficienti temporali si distribuisce nel dominio della frequenza.

## 10. Misure spettrali e di coerenza

Le misure di coerenza quantificano l'associazione lineare tra segnali nel dominio della frequenza. Sono utili per identificare frequenze condivise, ma non stabiliscono da sole una direzione causale.

### 10.1 Spettro, coerenza e coerenza immaginaria

| Stimatore | Formula | Riferimento e interpretazione |
|---|---|---|
| Matrice di densità spettrale | $\displaystyle \mathbf{S}(f)=\operatorname{E}\!\left[\mathbf{X}(f)\mathbf{X}(f)^{\mathrm H}\right]=\mathbf{H}(f)\boldsymbol{\Sigma}\mathbf{H}(f)^{\mathrm H}$ | Brillinger (2001). $S_{ii}(f)$ è lo spettro della variabile $i$; $S_{ij}(f)$ è il cross-spettro tra $i$ e $j$. |
| Coerenza (Coh) | $\displaystyle C_{ij}(f)=\frac{S_{ij}(f)}{\sqrt{S_{ii}(f)S_{jj}(f)}}$, con $\displaystyle 0\leq\lvert C_{ij}(f)\rvert^2\leq1$ | Brillinger (2001). È l'analogo in frequenza della correlazione normalizzata. $\lvert C_{ij}(f)\rvert^2$ è la magnitude-squared coherence; la fase può suggerire un ritardo, ma richiede cautela. |
| Coerenza immaginaria (iCoh) | $\displaystyle iCoh_{ij}(f)=\operatorname{Im}\!\left[C_{ij}(f)\right]$ | Nolte et al. (2004). Considera la parte immaginaria della coerenza e riduce il contributo degli accoppiamenti istantanei a fase nulla, come quelli spesso associati alla conduzione di volume. |

La coerenza è simmetrica in ampiezza e non identifica una sorgente e un bersaglio. Una fase diversa da zero può essere compatibile con un ritardo, ma non dimostra una relazione causale. La fase può dipendere anche dal filtraggio, dal riferimento, dal mixing e da sorgenti comuni.

La coerenza immaginaria attenua le componenti a fase nulla, ma non elimina ogni effetto del mixing. Può inoltre scartare relazioni fisiologiche realmente presenti a fase prossima a zero.

### 10.2 Coerenza condizionata e multipla

Per tenere conto delle altre variabili osservate si possono utilizzare la coerenza parziale e la coerenza multipla.

Definiamo la matrice spettrale inversa

$$
\widetilde{\mathbf{S}}(f)
=
\mathbf{S}(f)^{-1}.
$$

| Stimatore | Formula | Riferimento e interpretazione |
|---|---|---|
| Coerenza parziale (pCoh) | $\displaystyle P_{ij}(f)=-\frac{\widetilde{S}_{ij}(f)}{\sqrt{\widetilde{S}_{ii}(f)\widetilde{S}_{jj}(f)}}$, con $\displaystyle 0\leq\lvert P_{ij}(f)\rvert^2\leq1$ | Brillinger (2001). Misura la coerenza residua tra $i$ e $j$ dopo avere rimosso la parte spiegabile linearmente dalle altre variabili misurate. È una misura condizionata, ma non direzionale. |
| Coerenza multipla (mCoh) | $\displaystyle G_i(f)=\sqrt{1-\frac{\det\mathbf{S}(f)}{S_{ii}(f)M_{ii}(f)}}$ | Brillinger (2001). Quantità univariata che misura la coerenza complessiva della variabile $i$ con tutte le altre variabili osservate. $M_{ii}(f)$ è il minore di $\mathbf{S}(f)$ ottenuto eliminando la riga e la colonna $i$. |

Il condizionamento dipende dall'insieme di variabili presenti nel modello. Una coerenza parziale ridotta non esclude l'esistenza di un legame fisiologico; indica che l'associazione lineare osservata può essere spiegata, almeno in parte, dalle altre variabili misurate.

## 11. Partial Directed Coherence

La Partial Directed Coherence è costruita a partire dalla matrice $\mathbf{A}(f)$. Con la convenzione utilizzata nelle slide, l'indice $j$ identifica la sorgente e l'indice $i$ il bersaglio. Pertanto $\pi_{ij}(f)$ rappresenta la direzione

$$
j\rightarrow i.
$$

La PDC normalizza la connessione $A_{ij}(f)$ rispetto a tutte le connessioni uscenti dalla sorgente $j$. Risponde quindi alla domanda:

> Come si distribuisce l'influenza diretta della sorgente $j$ tra i possibili bersagli del sistema?

### 11.1 PDC e varianti

| Stimatore | Formula | Riferimento e interpretazione |
|---|---|---|
| PDC normalizzata | $\displaystyle \pi_{ij}(f)=\frac{A_{ij}(f)}{\sqrt{\sum_{k=1}^{M}\lvert A_{kj}(f)\rvert^2}}$, con $\displaystyle 0\leq\lvert\pi_{ij}(f)\rvert^2\leq1$ e $\displaystyle\sum_{i=1}^{M}\lvert\pi_{ij}(f)\rvert^2=1$ | Baccalá e Sameshima (2001). Misura complessa interpretabile come influenza condizionale da $j$ a $i$, normalizzata rispetto al flusso causale complessivo uscente da $j$. In genere si rappresenta $\lvert\pi_{ij}(f)\rvert^2$. |
| PDC generalizzata (gPDC) | $\displaystyle \overline{\pi}_{ij}(f)=\frac{A_{ij}(f)/\sigma_i}{\sqrt{\sum_{k=1}^{M}\lvert A_{kj}(f)\rvert^2/\sigma_k^2}}$, con $\displaystyle \sigma_i^2=\Sigma_{ii}$ | Baccalá e Sameshima (2007). Pesa i coefficienti mediante la deviazione standard delle innovazioni e riduce la dipendenza della PDC da scale o livelli di rumore diversi. |
| PDC rinormalizzata (rPDC) | $\displaystyle \lambda_{ij}(f)=\mathbf{Q}_{ij}(f)^{\mathsf T}\mathbf{V}_{ij}(f)^{-1}\mathbf{Q}_{ij}(f)$, con $\displaystyle \mathbf{Q}_{ij}(f)=\begin{bmatrix}\operatorname{Re}A_{ij}(f)\\ \operatorname{Im}A_{ij}(f)\end{bmatrix}$ e $\displaystyle \mathbf{V}_{ij}(f)=\operatorname{Cov}\!\left[\mathbf{Q}_{ij}(f)\right]$ | Schelter et al. (2009). Usa la covarianza della stima dei coefficienti VAR per ottenere una statistica indipendente dall'unità di misura e adatta all'inferenza. |

La normalizzazione della PDC avviene per colonna di $\mathbf{A}(f)$. A sorgente $j$ fissata, il numeratore seleziona la connessione $j\rightarrow i$, mentre il denominatore contiene tutte le connessioni che partono da $j$.

![Normalizzazione della PDC rispetto ai flussi uscenti dalla sorgente](../assets/figures/pdc_normalizzazione.png)

*Nella PDC la connessione $j\rightarrow i$ viene confrontata con tutte le connessioni uscenti dalla sorgente $j$.*

La PDC viene interpretata come misura diretta perché dipende da $A_{ij}(f)$. Se tutti i coefficienti $A_{k,ij}$ del VAR sono nulli, allora anche $A_{ij}(f)$ e la corrispondente PDC sono nulli. Il termine diretto resta però relativo alle variabili incluse nel modello: una variabile intermedia non osservata può produrre una connessione apparente tra i nodi misurati.

La PDC è una misura relativa. Se la sorgente $j$ distribuisce la propria influenza tra più bersagli, la quota attribuita a ciascuna connessione può diminuire anche senza una riduzione del coefficiente non normalizzato. Il valore non deve quindi essere interpretato come una quantità assoluta di potenza trasferita.

La gPDC è particolarmente utile quando i segnali o le innovazioni presentano scale molto diverse. La rPDC risponde invece a un'esigenza inferenziale: tiene conto dell'incertezza associata alla stima dei coefficienti.

## 12. Directed Transfer Function

La Directed Transfer Function è costruita a partire dalla matrice di trasferimento $\mathbf{H}(f)$. Anche in questo caso l'elemento con indici $ij$ rappresenta la direzione

$$
j\rightarrow i.
$$

La DTF normalizza $H_{ij}(f)$ rispetto a tutti i contributi che raggiungono il bersaglio $i$. Risponde quindi alla domanda:

> Quale quota dell'attività che raggiunge il bersaglio $i$ proviene dalla sorgente $j$?

Poiché $\mathbf{H}(f)=\mathbf{A}(f)^{-1}$ descrive la propagazione attraverso l'intero sistema, la DTF misura un flusso totale. Può quindi essere diversa da zero anche quando il percorso dalla sorgente al bersaglio attraversa una o più variabili intermedie.

### 12.1 DTF e varianti

| Stimatore | Formula | Riferimento e interpretazione |
|---|---|---|
| Directed Transfer Function (DTF) | $\displaystyle \gamma_{ij}(f)=\frac{H_{ij}(f)}{\sqrt{\sum_{k=1}^{M}\lvert H_{ik}(f)\rvert^2}}$, con $\displaystyle 0\leq\lvert\gamma_{ij}(f)\rvert^2\leq1$ | Kamiński e Blinowska (1991). Misura il flusso informativo totale da $j$ a $i$, normalizzato rispetto all'afflusso complessivo verso il bersaglio $i$. |
| Full-frequency DTF (ffDTF) | $\displaystyle \eta_{ij}^{2}(f)=\frac{\lvert H_{ij}(f)\rvert^2}{\sum_{\nu}\sum_{k=1}^{M}\lvert H_{ik}(\nu)\rvert^2}$ | Korzeniewska et al. (2003). Il denominatore integra su tutte le frequenze $\nu$ e facilita il confronto del flusso informativo tra frequenze diverse. |
| Direct DTF (dDTF) | $\displaystyle \delta_{ij}^{2}(f)=\eta_{ij}^{2}(f)\lvert P_{ij}(f)\rvert^2$ | Korzeniewska et al. (2003). Combina ffDTF e coerenza parziale; attenua i percorsi indiretti ed evidenzia la componente diretta del flusso direzionato. |

La normalizzazione della DTF avviene per riga di $\mathbf{H}(f)$. A bersaglio $i$ fissato, il numeratore seleziona il contributo $j\rightarrow i$, mentre il denominatore contiene tutti i contributi che raggiungono $i$.

![Normalizzazione della DTF rispetto agli afflussi verso il bersaglio](../assets/figures/dtf_normalizzazione.png)

*Nella DTF la connessione $j\rightarrow i$ viene confrontata con tutti i contributi che raggiungono il bersaglio $i$.*

Consideriamo il percorso

$$
1\rightarrow2\rightarrow3.
$$

Un'innovazione generata nel nodo 1 può propagarsi fino al nodo 3 attraverso il nodo 2. Di conseguenza, $H_{31}(f)$ e la DTF da 1 a 3 possono essere diversi da zero anche in assenza di un coefficiente VAR diretto da 1 a 3. Questo comportamento non rappresenta un errore: corrisponde alla definizione di trasferimento totale.

La ffDTF modifica la normalizzazione integrando il denominatore su tutte le frequenze. La dDTF moltiplica la ffDTF per la coerenza parziale e cerca così di attenuare i contributi dovuti ai percorsi indiretti.

### 12.2 Confronto tra PDC e DTF

| Aspetto | PDC | DTF |
|---|---|---|
| Matrice di partenza | $\mathbf{A}(f)$ | $\mathbf{H}(f)=\mathbf{A}(f)^{-1}$ |
| Elemento per $j\rightarrow i$ | $A_{ij}(f)$ | $H_{ij}(f)$ |
| Normalizzazione | Rispetto alle uscite della sorgente $j$ | Rispetto agli ingressi del bersaglio $i$ |
| Significato principale | Influenza diretta e condizionata rispetto alle variabili osservate | Trasferimento totale, comprendente anche percorsi indiretti |
| Domanda | Verso quali bersagli si distribuisce l'influenza della sorgente $j$? | Da quali sorgenti proviene l'attività che raggiunge il bersaglio $i$? |

PDC e DTF non rappresentano due stime intercambiabili della stessa quantità. Una connessione può risultare evidente nella DTF perché trasmette attività attraverso un percorso indiretto e risultare assente nella PDC, che riflette i coefficienti diretti del VAR.

Anche la normalizzazione influenza il confronto. Aggiungere o rimuovere una variabile dal modello modifica l'insieme delle uscite considerate dalla PDC e degli ingressi considerati dalla DTF. I valori normalizzati possono quindi cambiare anche quando il coefficiente relativo a una singola coppia rimane simile.

La scelta della misura deve seguire la domanda fisiologica:

- la PDC è adatta quando interessa la distribuzione delle influenze dirette in uscita da una sorgente;
- la DTF è adatta quando interessa ricostruire l'origine dell'attività che raggiunge un bersaglio, anche attraverso percorsi intermedi;
- la lettura congiunta di PDC e DTF aiuta a distinguere una connessione diretta da una propagazione indiretta.

Qualunque misura deve essere interpretata insieme agli spettri di potenza e alla significatività statistica. Un valore normalizzato elevato in una banda con potenza trascurabile può avere un'importanza fisiologica limitata.

## 13. Inferenza statistica per PDC e DTF

Le stime di PDC e DTF ottenute da un campione finito sono quasi sempre maggiori di zero, anche quando nel processo che ha generato i dati non esiste una reale interazione nella direzione considerata. Questo accade perché i coefficienti del modello VAR sono stimati con incertezza e tale incertezza si trasferisce alle misure di connettività.

L'inferenza statistica serve quindi a stabilire se il valore osservato di una misura sia sufficientemente grande rispetto alla variabilità attesa sotto un'opportuna ipotesi nulla. Indicando genericamente con

$$
\theta_{i\leftarrow j}(f)
$$

la misura di connettività dalla variabile sorgente $j$ alla variabile bersaglio $i$ alla frequenza $f$, il test è

$$
H_0:\theta_{i\leftarrow j}(f)=0
\qquad \text{contro} \qquad
H_1:\theta_{i\leftarrow j}(f)>0.
$$

Il significato di $H_0$ dipende però dalla misura utilizzata:

| Misura | Significato dell'ipotesi nulla |
|---|---|
| PDC | assenza di una connessione **diretta** da $j$ a $i$ |
| DTF | assenza di un trasferimento **complessivo**, diretto o indiretto, da $j$ a $i$ alla frequenza $f$ |

Questa distinzione è essenziale: una PDC non significativa e una DTF non significativa non esprimono esattamente la stessa affermazione sulla rete.

### 13.1 Distribuzioni asintotiche

Per campioni sufficientemente grandi, lo stimatore dei coefficienti VAR può essere approssimato mediante una distribuzione normale multivariata. Se $\mathbf a$ raccoglie tutti i coefficienti del modello e $\widehat{\mathbf a}$ è il corrispondente stimatore, allora

$$
\sqrt{N}\left(\widehat{\mathbf a}-\mathbf a\right)
\overset{a}{\sim}
\mathcal N\left(\mathbf 0,\mathbf V_a\right),
$$

dove $N$ è il numero effettivo di osservazioni e $\mathbf V_a$ è la matrice di covarianza asintotica dello stimatore.

Nel dominio della frequenza, per $i\neq j$, l'elemento fuori diagonale della matrice del modello è

$$
A_{ij}(f)
=
-\sum_{k=1}^{p} A_{ij}^{(k)}
e^{-\mathrm{i}2\pi f k\Delta t}.
$$

Si tratta di una combinazione lineare dei coefficienti VAR; di conseguenza, la sua parte reale e la sua parte immaginaria hanno una distribuzione congiuntamente gaussiana in senso asintotico. PDC e DTF sono però funzioni non lineari dei coefficienti, perché contengono moduli quadratici e normalizzazioni. Per esempio,

$$
\left|A_{ij}(f)\right|^2
=
\operatorname{Re}\!\left\{A_{ij}(f)\right\}^2
+
\operatorname{Im}\!\left\{A_{ij}(f)\right\}^2.
$$

Sotto $H_0$, la distribuzione asintotica della misura quadratica può essere espressa nella forma generale

$$
N\,\theta_{i\leftarrow j}(f)
\overset{a}{\sim}
\sum_r \lambda_r(f)\,\chi^2_{1,r},
$$

dove i pesi $\lambda_r(f)$ dipendono dal modello stimato e dalla misura considerata. Da questa distribuzione si può ricavare una soglia critica per ogni direzione e frequenza.

L'approssimazione asintotica è più affidabile quando:

- $N$ è grande rispetto al numero di parametri, che cresce come $M^2p$;
- il modello VAR è stabile e l'ordine $p$ è adeguato;
- i residui sono compatibili con le ipotesi del modello;
- il processo non è vicino al limite di instabilità.

Può invece essere poco accurata con campioni brevi, molti canali, ordine elevato, specificazione errata del modello o residui fortemente non gaussiani o correlati. Inoltre, la soglia ottenuta in questo modo è normalmente **puntuale**: controlla l'errore per una singola direzione e una singola frequenza, non per l'intera famiglia di confronti.

### 13.2 Inferenza mediante dati surrogati

Quando non si vuole affidare l'inferenza esclusivamente all'approssimazione asintotica, è possibile costruire empiricamente la distribuzione nulla mediante dati surrogati. L'idea generale è la seguente:

1. si definisce con precisione l'ipotesi nulla da verificare;
2. si generano $B$ dataset surrogati compatibili con tale ipotesi;
3. per ogni dataset si ristima il modello VAR;
4. dal modello ristimato si calcolano nuovamente PDC o DTF;
5. si confronta la misura osservata con la distribuzione delle misure surrogate.

Indicando con $\theta_{i\leftarrow j}^{*(b)}(f)$ la misura ottenuta dal surrogato $b$ e con $\theta_{i\leftarrow j}^{\mathrm{obs}}(f)$ quella osservata, un $p$-value empirico può essere calcolato come

$$
p
=
\frac{
1+\displaystyle\sum_{b=1}^{B}
\mathbb I\!\left[
\theta_{i\leftarrow j}^{*(b)}(f)
\geq
\theta_{i\leftarrow j}^{\mathrm{obs}}(f)
\right]
}{B+1}.
$$

Il termine $+1$ evita di assegnare un $p$-value nullo quando nessun surrogato supera il valore osservato. Ne segue che la risoluzione minima ottenibile è $1/(B+1)$: il numero di surrogati deve quindi essere scelto anche in funzione del livello di significatività desiderato.

#### 13.2.1 Randomizzazione delle fasi

La *phase randomization* costruisce surrogati che mantengono la dinamica lineare di ciascun segnale preso singolarmente, ma eliminano le relazioni coerenti tra segnali diversi.

Per ogni variabile $x_i(t)$ si calcola la trasformata di Fourier

$$
X_i(f)=|X_i(f)|e^{\mathrm{i}\phi_i(f)},
$$

si sostituisce la fase originale con una fase casuale indipendente

$$
\phi_i^*(f)\sim\mathcal U(0,2\pi),
$$

e si costruisce

$$
X_i^*(f)=|X_i(f)|e^{\mathrm{i}\phi_i^*(f)}.
$$

Le fasi devono rispettare la simmetria coniugata, così che la trasformata inversa produca un segnale reale. Le fasi vengono randomizzate indipendentemente per ciascun canale; su ogni dataset così ottenuto si ristimano il VAR e le misure di connettività.

| La randomizzazione delle fasi preserva | La randomizzazione delle fasi elimina |
|---|---|
| spettro di potenza di ciascuna serie | spettro incrociato tra serie |
| autocorrelazione univariata | relazioni di fase tra canali |
| media, se la componente continua non viene modificata, e varianza | dipendenze direzionali e non direzionali tra variabili |
| dinamica lineare univariata | struttura multivariata originale |

La procedura rappresenta quindi una **ipotesi nulla globale di assenza di dipendenza tra le variabili**. Non verifica selettivamente l'assenza della sola connessione $j\to i$ e può risultare troppo distruttiva quando si vuole testare un singolo arco conservando il resto della rete. La [Demo 51](../demo/demo_51_surrogati_fase.md) mostra operativamente la costruzione dei surrogati e il calcolo delle soglie.

#### 13.2.2 Causal shuffling

Il *causal shuffling* costruisce invece una distribuzione nulla modificando in modo selettivo la struttura causale del modello MVAR stimato. In termini generali, la procedura consiste nel:

1. stimare il modello MVAR sui dati originali;
2. trasformare il modello in modo da annullare la componente causale sottoposta a test;
3. simulare nuovi dati dal modello modificato;
4. ristimare il VAR e la connettività su ciascuna simulazione;
5. ripetere la procedura $B$ volte per ottenere la distribuzione nulla.

Le slide distinguono due formulazioni, coerenti con il diverso significato di PDC e DTF:

| Procedura | Ipotesi nulla | Interpretazione | Misura associata |
|---|---|---|---|
| CFTd | $H_0^{\mathrm{direct}}: A_{ij}^{(1)}=\cdots=A_{ij}^{(p)}=0$ | assenza della connessione diretta $j\to i$ a tutti i ritardi | PDC |
| CFTf | $H_0^{\mathrm{full}}: H_{ij}(f)=0$ | assenza del trasferimento complessivo $j\to i$ alla frequenza $f$ | DTF |

Rispetto alla randomizzazione delle fasi, questo approccio cerca di preservare le altre dinamiche e connessioni della rete. L'ipotesi nulla è quindi più mirata, ma la procedura dipende dal modello stimato: se il VAR di partenza rappresenta male i dati, anche i surrogati erediteranno tale limite (Faes et al., 2010).

### 13.3 Soglie puntuali e confronti multipli

In una rete con $M$ variabili si possono esaminare numerose direzioni e, per ciascuna, molte frequenze. Si esegue quindi una famiglia di test

$$
H_{0,ij}(f):\theta_{i\leftarrow j}(f)=0.
$$

Applicare separatamente una soglia al 5% a ogni punto aumenta la probabilità di ottenere almeno un falso positivo. La correzione deve riflettere la famiglia di risultati sulla quale si intende formulare la conclusione: tutte le frequenze di una direzione, tutte le direzioni di una banda oppure l'intera matrice frequenza-per-direzione.

Una soluzione particolarmente naturale con i surrogati è la statistica del massimo. Per il surrogato $b$ si calcola

$$
T_b^{\max}
=
\max_{i\neq j,\,f}
\theta_{i\leftarrow j}^{*(b)}(f),
$$

e si determina la soglia globale

$$
q_{1-\alpha}^{\max}
=
\operatorname{quantile}
\left(
\left\{T_1^{\max},\ldots,T_B^{\max}\right\};
1-\alpha
\right).
$$

La misura osservata è dichiarata significativa quando

$$
\theta_{i\leftarrow j}^{\mathrm{obs}}(f)
>
q_{1-\alpha}^{\max}.
$$

| Metodo | Errore controllato o principio | Osservazioni |
|---|---|---|
| Bonferroni | probabilità di almeno un falso positivo nella famiglia | semplice, ma spesso conservativo |
| FDR | proporzione attesa di falsi positivi tra i risultati dichiarati significativi | utile quando i confronti sono numerosi |
| Max-statistic | probabilità di almeno un falso positivo nella famiglia | si integra direttamente con i surrogati e tiene conto della distribuzione del massimo |
| Test basati su cluster | significatività di regioni contigue in frequenza o tempo-frequenza | sfruttano la continuità, ma richiedono di definire soglia iniziale e statistica di cluster |

Le bande di frequenza e le connessioni di interesse dovrebbero essere definite, quando possibile, prima di osservare i risultati. Ridurre la famiglia di test dopo avere ispezionato le mappe introduce una selezione dipendente dai dati e rende più difficile interpretare il livello di errore effettivamente controllato.

### 13.4 Interpretazione e presentazione dei risultati

Un risultato significativo indica che la misura osservata è incompatibile con l'ipotesi nulla scelta al livello di errore stabilito. Non dimostra, da solo, l'esistenza di un meccanismo fisiologico diretto: l'interpretazione resta subordinata alla qualità del modello VAR, alle variabili incluse, al riferimento dei segnali e alla possibile presenza di sorgenti comuni non osservate.

Analogamente, il mancato rifiuto di $H_0$ non dimostra che la connessione sia assente. Può riflettere una potenza statistica insufficiente, un campione breve, una stima variabile o una correzione multipla molto severa.

Per rendere l'analisi riproducibile è opportuno riportare almeno:

- la misura analizzata e la convenzione degli indici $i\leftarrow j$;
- l'ipotesi nulla effettivamente implementata;
- il metodo inferenziale, asintotico o basato su surrogati;
- il numero $B$ di surrogati e il relativo schema di generazione;
- la famiglia di confronti e la correzione applicata;
- il livello $\alpha$ e il tipo di soglia, puntuale o globale;
- l'ordine del VAR, il numero di campioni e gli esiti della validazione del modello.

In sintesi, una soglia statistica è interpretabile soltanto insieme alla procedura che l'ha generata. Cambiare il tipo di surrogato significa cambiare l'ipotesi nulla; cambiare la famiglia di confronti significa cambiare l'errore controllato.

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

Il modello VAR introdotto nelle sezioni precedenti assume che la relazione tra le variabili rimanga costante durante tutto l'intervallo analizzato. In altre parole, un unico insieme di matrici $\mathbf A_1,\ldots,\mathbf A_p$ descrive l'intera registrazione. Questa ipotesi può essere ragionevole per un segmento breve e sufficientemente omogeneo, ma diventa restrittiva quando il sistema attraversa condizioni diverse.

Nei segnali fisiologici la dinamica può cambiare nel tempo per diverse ragioni. Una risposta evocata può produrre una variazione transitoria, un soggetto può passare da uno stato di riposo a un compito, oppure una rete neurale può modificare le proprie interazioni durante le diverse fasi di uno stimolo. Possono cambiare la media e la varianza dei singoli segnali, ma possono cambiare anche le relazioni predittive tra le variabili. Questi due fenomeni non sono equivalenti: una variazione dell'ampiezza di un canale non implica necessariamente una variazione della connettività, così come una variazione dei coefficienti fuori diagonale può avvenire senza un grande cambiamento della potenza dei segnali.

Per rappresentare una dinamica che evolve nel tempo si introducono coefficienti VAR dipendenti da $t$:

$$
\mathbf{x}(t)
=
\sum_{r=1}^{p}\mathbf{A}_r(t)\mathbf{x}(t-r)
+
\mathbf{u}(t).
$$

La matrice $\mathbf A_r(t)$ descrive, al tempo $t$, come i campioni osservati al ritardo $r$ contribuiscono al valore corrente. Gli elementi diagonali rappresentano la dipendenza di ciascuna variabile dal proprio passato, mentre gli elementi fuori diagonale descrivono i contributi ritardati tra variabili diverse. Se questi ultimi cambiano, cambia anche la struttura della rete descritta dal modello.

La dipendenza dal tempo non significa però che ogni coefficiente possa essere stimato indipendentemente a ogni campione. Un singolo campione non contiene informazione sufficiente per identificare tutti i parametri di un modello multivariato. Ogni metodo tempo-variante deve quindi raccogliere informazione da un intorno temporale oppure trasferirla dalle stime precedenti. In pratica si introduce sempre un'ipotesi di regolarità: i coefficienti devono rimanere approssimativamente costanti per un breve intervallo oppure devono evolvere in modo abbastanza graduale da poter essere seguiti da un algoritmo adattivo.

L'obiettivo non è più ottenere una sola rete media per l'intera registrazione, ma ricostruire una successione di modelli locali. Da ciascun modello si possono calcolare coefficienti, causalità di Granger, PDC o DTF dipendenti dal tempo. Per esempio, una misura spettrale diventa una funzione di tempo e frequenza,

$$
\theta_{i\leftarrow j}(t,f),
$$

che descrive come la relazione dalla sorgente $j$ al bersaglio $i$ varia nel corso dell'esperimento e alle diverse frequenze.

Questa maggiore flessibilità introduce un problema fondamentale. Una variazione osservata nella connettività stimata può riflettere un cambiamento reale del processo, ma può anche dipendere dal rumore, dal numero limitato di osservazioni o dalla rapidità con cui l'algoritmo è stato configurato per adattarsi. Una mappa molto variabile non costituisce quindi, da sola, una prova di dinamica fisiologica. Occorre verificare che il risultato sia stabile rispetto alle scelte del metodo e che i modelli locali rimangano adeguati ai dati.

I principali approcci sono:

- le **finestre mobili**, che stimano un modello separato su intervalli successivi;
- i metodi ricorsivi, come il **Recursive Least Squares**, che aggiornano i coefficienti quando arriva un nuovo campione;
- i modelli stato-spazio, nei quali i coefficienti costituiscono uno stato nascosto stimato, per esempio, mediante il **filtro di Kalman**.

Tutti introducono un compromesso tra rapidità di adattamento e affidabilità della stima. Se il metodo usa molta informazione passata, le stime risultano più regolari ma seguono con ritardo i cambiamenti. Se usa soprattutto dati recenti, localizza meglio le transizioni ma diventa più sensibile alla variabilità campionaria.

### 15.1 Finestre mobili

L'approccio più diretto consiste nel suddividere la registrazione in finestre temporali e stimare un modello VAR separato in ciascuna finestra. Se la finestra contiene $W$ campioni e viene spostata di $Q$ campioni a ogni iterazione, il numero di modelli stimati su una registrazione di $T$ campioni è

$$
N_w
=
\left\lfloor\frac{T-W}{Q}\right\rfloor+1.
$$

Per la finestra centrata intorno al tempo $t_\ell$ si ottengono le matrici locali

$$
\widehat{\mathbf A}_1(t_\ell),\ldots,
\widehat{\mathbf A}_p(t_\ell).
$$

Il valore attribuito a $t_\ell$ riassume quindi l'intero intervallo coperto dalla finestra. Non rappresenta una stima istantanea nel senso stretto del termine. Se una connessione cambia bruscamente all'interno della finestra, il modello combina campioni precedenti e successivi alla transizione e restituisce una dinamica media. Di conseguenza, un cambiamento realmente brusco può apparire più graduale nelle mappe stimate.

#### Stazionarietà locale

L'assunzione centrale è che, all'interno di ogni finestra, il processo sia approssimativamente stazionario e possa essere descritto da coefficienti costanti. La registrazione completa può dunque essere non stazionaria, purché le sue proprietà cambino abbastanza lentamente rispetto alla durata della finestra.

La stazionarietà locale rimane un'approssimazione. Una finestra che attraversa l'inizio di uno stimolo o una transizione di stato può contenere due regimi differenti. In quel caso il VAR stimato non descrive bene nessuno dei due regimi, ma una loro combinazione. È quindi utile mettere in relazione la posizione delle finestre con la struttura del protocollo sperimentale e con l'ispezione dei segnali.

#### Scelta della lunghezza della finestra

La lunghezza $W$ determina il principale compromesso del metodo:

| Scelta di $W$ | Vantaggio | Limite |
|---|---|---|
| finestra lunga | più osservazioni e coefficienti generalmente più stabili | le transizioni vengono mediate e la risoluzione temporale diminuisce |
| finestra corta | migliore localizzazione temporale | maggiore varianza e rischio di stimare un modello poco affidabile |

La finestra deve essere abbastanza lunga da contenere informazioni sulla dinamica che si vuole studiare. Per analisi nel dominio della frequenza deve inoltre comprendere un numero adeguato di cicli delle oscillazioni di interesse. Una finestra molto breve rispetto al periodo delle frequenze più basse non permette di caratterizzarle in modo attendibile.

La numerosità disponibile per la stima non coincide esattamente con $W$. Per un VAR di ordine $p$, i primi $p$ campioni servono a costruire i regressori, quindi ogni finestra fornisce circa $W-p$ osservazioni utilizzabili. All'aumentare del numero di variabili $M$ e dell'ordine $p$ cresce anche il numero di coefficienti da stimare. Una finestra accettabile per un modello bivariato di basso ordine può quindi essere insufficiente per una rete più grande.

#### Passo di scorrimento e sovrapposizione

Il passo $Q$ stabilisce ogni quanti campioni viene prodotta una nuova stima. Quando $Q<W$, le finestre si sovrappongono e condividono parte dei dati. Un passo piccolo genera una mappa più densamente campionata e visivamente più continua, ma non aumenta la quantità di informazione indipendente contenuta nella registrazione.

È importante distinguere:

- la **durata della finestra** $W\Delta t$, che determina quanto esteso è l'intervallo temporale riassunto da ciascun modello;
- il **passo** $Q\Delta t$, che determina la distanza temporale tra due stime consecutive.

Ridurre $Q$ non compensa una finestra troppo lunga. Per esempio, stimare un modello ogni 10 ms con finestre di 1 s produce valori ogni 10 ms, ma ciascun valore continua a riassumere circa 1 s di dati. La risoluzione nominale della mappa è quindi diversa dalla capacità effettiva di distinguere due cambiamenti ravvicinati.

Poiché finestre adiacenti condividono molti campioni, anche le stime sono fortemente dipendenti. Una sequenza regolare non contiene necessariamente molte conferme indipendenti dello stesso effetto. Questa dipendenza deve essere considerata soprattutto nell'inferenza sulle mappe tempo-frequenza.

#### Uso di più prove sperimentali

Quando sono disponibili più trial, si possono riunire, per ogni posizione temporale, i dati provenienti dalla stessa finestra di tutte le prove. Con $R$ trial, il numero complessivo di osservazioni utilizzabili diventa approssimativamente $R(W-p)$. Questa strategia permette di mantenere una finestra relativamente breve senza ridurre eccessivamente la precisione della stima.

Il pooling richiede però che i trial rappresentino realizzazioni ripetute dello stesso processo locale e che gli eventi siano allineati in modo adeguato. Jitter temporale, adattamento nel corso dell'esperimento o differenze sistematiche tra prove possono attenuare una variazione reale oppure generare una dinamica apparente. Prima di aggregare i trial è quindi necessario valutarne la comparabilità.

#### Scelte che rendono confrontabili le finestre

Per interpretare le differenze tra istanti come variazioni del processo conviene mantenere coerente la procedura di stima:

- usare lo stesso insieme di variabili e lo stesso ordine $p$ in tutte le finestre;
- applicare criteri di preprocessing comparabili;
- controllare localmente stabilità del VAR e bianchezza dei residui;
- verificare che i risultati non cambino sostanzialmente per piccole variazioni di $W$ e $Q$;
- segnalare le finestre nelle quali il rapporto tra osservazioni e parametri rende la stima poco affidabile.

Scegliere separatamente l'ordine del modello in ogni finestra può sembrare più flessibile, ma rende meno immediato il confronto temporale: una variazione della misura potrebbe derivare dal cambiamento di $p$ anziché dalla dinamica fisiologica. In genere è preferibile selezionare un ordine comune mediante un'analisi preliminare e valutarne poi l'adeguatezza nelle singole finestre.

```{figure} ../assets/figures/lunghezza_finestra.png
:alt: Compromesso tra finestra temporale corta e lunga
:width: 82%
:name: fig-lunghezza-finestra

La lunghezza della finestra determina il compromesso tra localizzazione temporale e affidabilità della stima. Il passo di scorrimento determina invece quanto fittamente vengono campionate le stime nel tempo. Schema tratto dalla presentazione del corso.
```

In sintesi, le finestre scorrevoli sono semplici da applicare e rendono esplicito quali dati contribuiscono a ciascuna stima. Il loro limite principale è che la scala temporale dell'analisi viene fissata a priori attraverso $W$. Le mappe ottenute devono pertanto essere interpretate come una sequenza di descrizioni locali mediate sulla finestra, non come misure istantanee della connettività.

La [Demo 52](../demo/demo_52_tv_var.md) applica le finestre mobili a una singola realizzazione multivariata e confronta il risultato con RLS e filtro di Kalman, usando gli stessi segnali e lo stesso ordine.

### 15.2 Recursive Least Squares

Il **Recursive Least Squares** (RLS) stima i coefficienti del VAR aggiornandoli ogni volta che diventa disponibile un nuovo campione. A differenza delle finestre mobili, non ristima il modello da zero su un blocco di dati: conserva una sintesi delle osservazioni precedenti e la modifica ricorsivamente.

Scriviamo il TV-VAR nella forma di regressione

$$
\mathbf{x}_t=\mathbf{B}_t\mathbf{z}_t+\mathbf{u}_t,
$$

dove

$$
\mathbf{B}_t=
\begin{bmatrix}
\mathbf{A}_1(t)&\cdots&\mathbf{A}_p(t)
\end{bmatrix}
\in\mathbb{R}^{M\times Mp}
$$

raccoglie tutti i coefficienti VAR e

$$
\mathbf{z}_t=
\begin{bmatrix}
\mathbf{x}_{t-1}^{\mathsf T}&
\cdots&
\mathbf{x}_{t-p}^{\mathsf T}
\end{bmatrix}^{\mathsf T}
\in\mathbb{R}^{Mp}
$$

contiene i campioni passati usati come regressori.

#### Minimi quadrati con forgetting factor esponenziale

Al tempo $t$, l'RLS cerca la matrice $\mathbf B$ che minimizza

$$
J_t(\mathbf B)
=
\sum_{s=p+1}^{t}
\lambda^{,t-s}
\left\|\mathbf{x}_s-\mathbf B\mathbf z_s\right\|_2^2,
\qquad 0<\lambda\leq1.
$$

Il peso di un'osservazione diminuisce geometricamente con la sua distanza dal tempo corrente. Il parametro $\lambda$ prende il nome di **forgetting factor**:

- con $\lambda=1$ tutte le osservazioni mantengono lo stesso peso e l'algoritmo ricorsivo converge alla stima ordinaria sui dati accumulati;
- con $\lambda<1$ le osservazioni recenti ricevono più peso e la stima può seguire coefficienti variabili nel tempo;
- diminuendo $\lambda$ aumenta la rapidità di adattamento, ma aumenta anche la variabilità della stima.

Una misura indicativa della memoria dell'algoritmo è

$$
N_{\mathrm{eff}}\approx\frac{1}{1-\lambda}.
$$

Questa quantità non rappresenta una finestra rigida: nessun campione viene eliminato bruscamente. Indica piuttosto l'ordine di grandezza del numero di osservazioni che contribuiscono maggiormente alla stima. Per esempio, valori di $\lambda$ molto vicini a uno producono una memoria lunga e un adattamento lento.

#### Aggiornamento ricorsivo

Prima dell'aggiornamento si calcola l'errore di previsione

$$
\mathbf e_t
=
\mathbf x_t-\widehat{\mathbf B}_{t-1}\mathbf z_t.
$$

Il vettore di guadagno è

$$
\mathbf{k}_t=
\frac{
\mathbf{P}_{t-1}\mathbf{z}_t
}{
\lambda+\mathbf{z}_t^{\mathsf T}
\mathbf{P}_{t-1}\mathbf{z}_t
},
$$

e i coefficienti vengono corretti secondo

$$
\widehat{\mathbf{B}}_t
=
\widehat{\mathbf{B}}_{t-1}
+
\mathbf{e}_t\mathbf{k}_t^{\mathsf T}.
$$

Infine si aggiorna la matrice

$$
\mathbf{P}_t
=
\frac{1}{\lambda}
\left(
\mathbf{P}_{t-1}
-
\mathbf{k}_t\mathbf{z}_t^{\mathsf T}
\mathbf{P}_{t-1}
\right).
$$

La matrice $\mathbf P_t$ può essere interpretata come una matrice di informazione inversa dei regressori. Se il nuovo vettore $\mathbf z_t$ contiene una configurazione già osservata molte volte, modifica poco la stima. Se porta informazione nuova rispetto alla storia recente, il guadagno aumenta e la correzione può essere maggiore.

L'aggiornamento $\mathbf e_t\mathbf k_t^{\mathsf T}$ ha dimensione $M\times Mp$, come $\mathbf B_t$. Lo stesso vettore di regressori viene usato per tutte le equazioni del VAR, mentre ciascuna componente di $\mathbf e_t$ corregge l'equazione della corrispondente variabile bersaglio.

#### Inizializzazione e interpretazione

L'algoritmo richiede una stima iniziale $\widehat{\mathbf B}_0$ e una matrice $\mathbf P_0$. I coefficienti possono essere inizializzati a zero oppure mediante un VAR stimato su un intervallo iniziale. Una scelta comune è

$$
\mathbf P_0=\delta^{-1}\mathbf I,
$$

dove un valore piccolo di $\delta$ esprime scarsa fiducia nella stima iniziale e permette correzioni più ampie. Le prime stime dipendono maggiormente dall'inizializzazione; è quindi opportuno considerare un periodo iniziale di assestamento e non interpretarlo come se avesse la stessa affidabilità del resto della registrazione.

L'RLS realizza un compromesso analogo a quello delle finestre mobili, ma sostituisce il bordo netto della finestra con una pesatura esponenziale. Una finestra assegna lo stesso peso a tutti i campioni interni e peso nullo a quelli esterni; l'RLS riduce progressivamente il peso delle osservazioni più lontane.

#### Limiti dell'RLS

Un unico valore di $\lambda$ impone la stessa scala di adattamento per tutta la registrazione e, nella formulazione standard, per tutti i coefficienti. Il metodo può quindi seguire male un processo che alterna lunghi periodi stabili e cambiamenti molto rapidi. Un valore troppo basso può far inseguire rumore e artefatti, mentre un valore troppo alto può attenuare o ritardare una transizione reale.

L'aggiornamento non garantisce inoltre che il VAR stimato sia stabile a ogni istante. Prima di calcolare misure spettrali occorre controllare la stabilità locale e osservare l'errore di previsione. Salti isolati dei coefficienti, soprattutto in corrispondenza di artefatti, non devono essere interpretati automaticamente come cambiamenti di connettività.

Nella [Demo 52](../demo/demo_52_tv_var.md), i coefficienti noti permettono di osservare il compromesso fra variabilità delle stime e ritardo nell'inseguimento delle transizioni.

## 16. Filtro di Kalman per il TV-VAR

Il filtro di Kalman fornisce una formulazione probabilistica della stima adattiva. I coefficienti variabili nel tempo vengono trattati come uno **stato nascosto**: non sono osservati direttamente, ma vengono inferiti dai segnali e dai loro valori passati. Il modello specifica sia come ci si aspetta che i coefficienti evolvano, sia come essi generano le osservazioni.

Il risultato, a ogni istante, comprende due elementi:

- una stima dei coefficienti;
- una matrice che descrive l'incertezza associata a tale stima.

Questa distinzione è importante. Due traiettorie dei coefficienti possono apparire simili, ma avere affidabilità molto diversa se una è sostenuta da regressori informativi e l'altra deriva da un intervallo rumoroso o scarsamente eccitato.

### 16.1 Modello stato-spazio generale

Un modello stato-spazio lineare è descritto da un'equazione di stato e da un'equazione di osservazione:

$$
\boldsymbol{\theta}_t
=
\mathbf F_t\boldsymbol{\theta}_{t-1}
+
\mathbf w_t,
\qquad
\mathbf w_t\sim
\mathcal N(\mathbf 0,\mathbf Q_{\theta}),
$$

$$
\mathbf y_t
=
\mathbf H_t\boldsymbol{\theta}_t
+
\mathbf v_t,
\qquad
\mathbf v_t\sim
\mathcal N(\mathbf 0,\mathbf R_y).
$$

L'equazione di stato descrive l'evoluzione dello stato nascosto $\boldsymbol\theta_t$. La matrice $\mathbf F_t$ rappresenta la dinamica attesa, mentre $\mathbf Q_{\theta}$ quantifica quanto lo stato può discostarsi da tale dinamica.

L'equazione di osservazione collega lo stato al dato misurato. La matrice $\mathbf H_t$ traduce i parametri nascosti nello spazio delle osservazioni e $\mathbf R_y$ descrive la variabilità che il modello di osservazione non spiega.

Si assume normalmente che $\mathbf w_t$ e $\mathbf v_t$ abbiano media nulla, siano serialmente incorrelati e indipendenti tra loro. La gaussianità rende esatta la descrizione mediante media e covarianza nel modello lineare. Anche quando essa è soltanto approssimata, il filtro conserva l'interpretazione di stimatore lineare ricorsivo.

### 16.2 Predizione e correzione

La notazione $t|s$ indica una quantità relativa al tempo $t$ stimata usando i dati disponibili fino al tempo $s$. Il filtro alterna una fase di predizione, basata sul modello di stato, e una fase di correzione, basata sulla nuova osservazione.

#### Predizione

La stima precedente viene propagata al tempo corrente:

$$
\widehat{\boldsymbol{\theta}}_{t|t-1}
=
\mathbf F_t
\widehat{\boldsymbol{\theta}}_{t-1|t-1}.
$$

Anche la covarianza dell'errore viene propagata:

$$
\mathbf P_{t|t-1}
=
\mathbf F_t
\mathbf P_{t-1|t-1}
\mathbf F_t^{\mathsf T}
+
\mathbf Q_{\theta}.
$$

Il termine $\mathbf Q_{\theta}$ aumenta l'incertezza perché ammette che lo stato possa essere cambiato tra $t-1$ e $t$.

#### Innovazione

La previsione dell'osservazione è

$$
\widehat{\mathbf y}_{t|t-1}
=
\mathbf H_t
\widehat{\boldsymbol{\theta}}_{t|t-1}.
$$

La differenza tra dato osservato e dato previsto prende il nome di innovazione del filtro:

$$
\boldsymbol{\nu}_t
=
\mathbf y_t-
\mathbf H_t\widehat{\boldsymbol{\theta}}_{t|t-1}.
$$

La sua covarianza prevista è

$$
\mathbf S_t
=
\mathbf H_t
\mathbf P_{t|t-1}
\mathbf H_t^{\mathsf T}
+
\mathbf R_y.
$$

Un'innovazione grande non implica necessariamente un cambiamento dello stato. Va confrontata con $\mathbf S_t$: uno scarto può essere rilevante quando l'incertezza prevista è piccola e del tutto plausibile quando osservazione e stato sono molto incerti.

#### Guadagno e correzione

Il guadagno di Kalman è

$$
\mathbf K_t
=
\mathbf P_{t|t-1}
\mathbf H_t^{\mathsf T}
\mathbf S_t^{-1}.
$$

La stima corretta diventa

$$
\widehat{\boldsymbol{\theta}}_{t|t}
=
\widehat{\boldsymbol{\theta}}_{t|t-1}
+
\mathbf K_t\boldsymbol{\nu}_t,
$$

e la sua covarianza viene aggiornata come

$$
\mathbf P_{t|t}
=
(\mathbf I-\mathbf K_t\mathbf H_t)
\mathbf P_{t|t-1}.
$$

Il guadagno bilancia due fonti di informazione. Se la previsione dello stato è incerta rispetto al rumore osservativo, il nuovo dato produce una correzione ampia. Se l'osservazione è molto rumorosa, il filtro si affida maggiormente alla traiettoria prevista.

| Situazione | Effetto prevalente |
|---|---|
| $\mathbf P_{t\|t-1}$ grande | maggiore disponibilità a correggere lo stato |
| $\mathbf R_y$ grande | minore peso assegnato alla nuova osservazione |
| $\mathbf Q_{\theta}$ grande | crescita più rapida dell'incertezza e maggiore adattabilità |
| osservazione informativa | riduzione della covarianza dopo la correzione |

Il filtro usa soltanto dati presenti e passati. Un **Kalman smoother**, applicato dopo avere acquisito l'intera registrazione, usa anche osservazioni future per rivedere le stime precedenti. Le mappe filtrate e quelle smussate non hanno quindi la stessa interpretazione temporale e il metodo utilizzato deve essere dichiarato.

### 16.3 Il TV-VAR come modello stato-spazio

Nel TV-VAR, lo stato nascosto raccoglie i coefficienti delle matrici autoregressive:

$$
\boldsymbol{\theta}_t
=
\operatorname{vec}
\begin{bmatrix}
\mathbf A_1(t)&\cdots&\mathbf A_p(t)
\end{bmatrix},
\qquad
\dim(\boldsymbol{\theta}_t)=M^2p.
$$

Il vettore dei regressori è

$$
\mathbf z_t
=
\begin{bmatrix}
\mathbf x_{t-1}^{\mathsf T}&
\cdots&
\mathbf x_{t-p}^{\mathsf T}
\end{bmatrix}^{\mathsf T}
\in\mathbb R^{Mp}.
$$

Usando il prodotto di Kronecker, l'equazione VAR può essere riscritta come

$$
\mathbf x_t
=
\underbrace{
(\mathbf z_t^{\mathsf T}\otimes\mathbf I_M)
}_{\mathbf H_t}
\boldsymbol{\theta}_t
+
\mathbf u_t,
$$

con

$$
\mathbf H_t\in
\mathbb R^{M\times M^2p}.
$$

La matrice di osservazione non è costante: cambia a ogni campione perché viene costruita con i valori passati dei segnali. In questa formulazione, $\mathbf u_t$ svolge il ruolo di rumore dell'equazione di osservazione e la sua covarianza è

$$
\mathbf R_u=\operatorname{Cov}(\mathbf u_t).
$$

Un modello semplice per l'evoluzione dei coefficienti è il random walk:

$$
\boldsymbol{\theta}_t
=
\boldsymbol{\theta}_{t-1}
+
\mathbf w_t,
$$

per cui $\mathbf F_t=\mathbf I$. Il random walk non afferma che i coefficienti fisiologici compiano realmente un moto casuale. Esprime una forma di regolarizzazione temporale: la migliore previsione iniziale è che i coefficienti persistano, ma $\mathbf Q_{\theta}$ permette loro di cambiare quando i dati forniscono evidenza sufficiente.

È utile distinguere due quantità chiamate spesso “innovazione”:

- l'innovazione VAR $\mathbf u_t$ rappresenta la parte di $\mathbf x_t$ non prevedibile mediante i ritardi nel modello vero;
- l'innovazione del Kalman $\boldsymbol\nu_t$ è lo scarto calcolato usando la stima predetta dei coefficienti e comprende anche la loro incertezza.

Le due quantità sono collegate, ma non sono intercambiabili. In particolare, la covarianza $\mathbf S_t$ contiene sia l'incertezza sui coefficienti sia $\mathbf R_u$.

La dimensione dello stato cresce come $M^2p$ e quella di $\mathbf P_t$ cresce con il quadrato della dimensione dello stato. L'aumento del numero di canali o dell'ordine produce quindi un costo numerico e statistico molto rapido. Formulazioni strutturate o generalizzate del filtro possono ridurre tale costo, ma richiedono ipotesi sulla struttura delle covarianze e sui coefficienti (Milde et al., 2010).

### 16.4 Inizializzazione e ruolo delle covarianze

Il comportamento del filtro dipende in modo sostanziale da quattro elementi:

| Quantità | Significato | Effetto di valori maggiori |
|---|---|---|
| $\widehat{\boldsymbol\theta}_{0\|0}$ | stima iniziale dei coefficienti | dipende dai valori assegnati, non da una scala unica |
| $\mathbf P_{0\|0}$ | incertezza iniziale sui coefficienti | maggiore correzione nelle prime osservazioni |
| $\mathbf Q_{\theta}$ | variabilità ammessa nell'evoluzione dei coefficienti | adattamento più rapido e traiettorie meno regolari |
| $\mathbf R_u$ | variabilità delle innovazioni VAR | minore peso del singolo campione osservato |

I coefficienti iniziali possono essere ottenuti da un VAR stazionario stimato su un intervallo iniziale oppure posti a zero quando non è disponibile una stima preliminare. La matrice $\mathbf P_{0|0}$ deve riflettere la fiducia in tale scelta. Un'inizializzazione incerta genera un transitorio durante il quale le stime possono cambiare rapidamente; questo intervallo va identificato e trattato con cautela.

Se si assume

$$
\mathbf Q_{\theta}=q\mathbf I,
$$

il parametro $q$ impone la stessa variabilità a tutti i coefficienti. Un valore piccolo produce traiettorie regolari ma può attenuare o ritardare cambiamenti rapidi. Un valore grande permette un adattamento veloce, ma aumenta il rischio che il filtro interpreti rumore e artefatti come variazioni dei coefficienti.

Una matrice $\mathbf Q_{\theta}$ non diagonale o strutturata può rappresentare velocità diverse e variazioni coordinate tra gruppi di coefficienti. Questa flessibilità richiede però più parametri e scelte più difficili da giustificare. La struttura adottata deve quindi essere riportata esplicitamente.

La matrice $\mathbf R_u$ descrive la variabilità non spiegata dal TV-VAR. Gli elementi fuori diagonale rappresentano variazioni contemporanee comuni alle equazioni, non connessioni ritardate. $\mathbf R_u$ può essere fissata usando un modello preliminare oppure aggiornata nel tempo, ma un aggiornamento troppo rapido rende difficile distinguere un cambiamento dei coefficienti da un cambiamento del rumore.

Il rapporto tra $\mathbf Q_{\theta}$ e $\mathbf R_u$ determina la scala temporale effettiva dell'adattamento:

- aumentando $\mathbf Q_{\theta}$ rispetto a $\mathbf R_u$, il filtro attribuisce più facilmente gli errori di previsione a un cambiamento dei coefficienti;
- aumentando $\mathbf R_u$ rispetto a $\mathbf Q_{\theta}$, il filtro li attribuisce più facilmente alla variabilità delle osservazioni.

Non esiste quindi una risoluzione temporale determinata soltanto dalla frequenza di campionamento. Anche se il filtro produce una stima a ogni campione, due cambiamenti ravvicinati possono essere indistinguibili quando la memoria implicita del filtro è lunga.

#### Dati con più trial

Con più prove sperimentali si può assumere che, allo stesso tempo relativo all'evento, i trial condividano le stesse matrici $\mathbf A_r(t)$. Le osservazioni dei diversi trial forniscono allora più informazione sullo stesso stato. Possono essere impilate in un'unica equazione di osservazione oppure elaborate con aggiornamenti equivalenti, purché sia specificata la struttura della covarianza tra trial.

Questa strategia aumenta l'informazione disponibile, ma richiede trial comparabili e correttamente allineati. Se la latenza della risposta varia tra prove, la stima può smussare la dinamica reale. Anche la standardizzazione deve essere descritta indicando lungo quale dimensione vengono calcolate media e deviazione standard.

#### Scelta pratica dei parametri

I parametri adattivi non dovrebbero essere scelti soltanto perché producono mappe visivamente nette. È preferibile combinarne la valutazione mediante:

- simulazioni con traiettorie dei coefficienti note;
- errore di previsione e analisi delle innovazioni;
- controllo della stabilità locale del VAR;
- analisi di sensibilità a $\mathbf Q_{\theta}$, $\mathbf R_u$ e inizializzazione;
- verifica della riproducibilità tra trial, soggetti o intervalli comparabili.

Una connessione che compare soltanto per una scelta molto specifica dei parametri adattivi è meno robusta di un risultato che persiste in un intervallo plausibile di configurazioni.

### 16.5 Dalla stima dei coefficienti alle mappe tempo-frequenza

Per ogni istante si ricostruiscono le matrici

$$
\widehat{\mathbf A}_1(t),\ldots,
\widehat{\mathbf A}_p(t)
$$

e la rappresentazione del modello nel dominio della frequenza:

$$
\widehat{\mathbf A}(f,t)
=
\mathbf I-
\sum_{r=1}^{p}
\widehat{\mathbf A}_r(t)
e^{-\mathrm{i}2\pi fr\Delta t}.
$$

La funzione di trasferimento locale è

$$
\widehat{\mathbf H}(f,t)
=
\widehat{\mathbf A}(f,t)^{-1}.
$$

Applicando a queste matrici le definizioni introdotte nelle sezioni precedenti si ottengono causalità di Granger tempo-variante, PDC e DTF tempo-frequenza. La notazione generale

$$
\theta_{i\leftarrow j}(t,f)
$$

indica la misura dalla sorgente $j$ al bersaglio $i$, valutata al tempo $t$ e alla frequenza $f$.

Prima dell'inversione di $\widehat{\mathbf A}(f,t)$ bisogna verificare la stabilità del modello locale. Un VAR instabile o molto vicino al limite di stabilità può generare picchi spettrali e valori di connettività estremi che dipendono dalla stima, non da una dinamica fisiologica plausibile.

La presenza di una stima per ogni campione non implica una risoluzione fisiologica campione per campione. Ogni valore dipende dalla memoria del filtro, dall'inizializzazione, dalle covarianze e dalla quantità di informazione contenuta nei regressori. Le mappe devono quindi essere descritte come stime adattive localizzate nel tempo, specificando i parametri che ne determinano la scala temporale effettiva.

È inoltre utile esaminare separatamente:

1. i coefficienti $\widehat{\mathbf A}_r(t)$;
2. la potenza o gli spettri locali dei segnali;
3. le misure normalizzate di connettività.

PDC e DTF possono cambiare anche per effetto dei rispettivi denominatori. Una diminuzione della PDC $j\to i$, per esempio, non implica necessariamente che il corrispondente coefficiente sia diminuito: potrebbero essere aumentate altre uscite dalla sorgente che partecipano alla normalizzazione. Affiancare coefficienti, spettri e connettività aiuta a distinguere questi casi.

### 16.6 Demo 52: simulazione TV-VAR e confronto dei tre metodi

La [Demo 52 — TV-VAR con finestre mobili, RLS e Kalman](../demo/demo_52_tv_var.md) riunisce i tre approcci in un unico esperimento. Si genera una registrazione di tre variabili, organizzata come matrice $3\times24000$, corrispondente a 120 s a 200 Hz. Non sono presenti trial ripetuti. I coefficienti noti consentono di distinguere i cambiamenti del processo dalla variabilità introdotta dalla stima.

Il modello è un VAR(2) con struttura a catena:

$$
x_1\longrightarrow x_2\longrightarrow x_3.
$$

Ogni variabile possiede una dinamica oscillatoria. La connessione $1\to2$ aumenta intorno a 30 s e torna al livello iniziale intorno a 78 s; la connessione $2\to3$ diminuisce intorno a 62 s. Le transizioni sono graduali. La frequenza propria del primo oscillatore passa inoltre da 10 a 15 Hz, con transizione centrata intorno a 55 s. Possiamo quindi osservare sia cambiamenti dell'intensità delle interazioni sia uno spostamento della loro distribuzione in frequenza.

La connessione diretta $1\to3$ rimane assente. Il segnale del primo nodo può propagarsi al terzo attraverso il secondo, ma la PDC teorica $1\to3$ è nulla perché descrive un arco diretto.

#### Stessi dati, tre strategie di adattamento

I tre metodi usano gli stessi segnali e lo stesso ordine $p=2$, fissato al valore vero per isolare il confronto degli stimatori.

| Metodo | Configurazione | Interpretazione |
|---|---|---|
| Finestre mobili | $W=1200$ campioni, pari a 6 s; passo $Q=20$ campioni, pari a 100 ms | ogni VAR riassume la finestra precedente |
| RLS | $\lambda=1-2/W\simeq0.998333$ | memoria esponenziale con scala indicativa di circa 3 s |
| Kalman | random walk con $\mathbf Q_\theta=2\cdot10^{-7}\mathbf I$ | adattamento controllato dalle covarianze di processo e delle innovazioni |

Le finestre vengono assegnate al loro **estremo destro**. Anche RLS e Kalman usano soltanto dati presenti e passati. Le traiettorie sono confrontate senza riallineamenti a posteriori: il ritardo nell'inseguimento delle transizioni rimane visibile.

RLS e Kalman condividono un'inizializzazione dei coefficienti ottenuta da un VAR sui primi 8 s. La covarianza $\mathbf R_u$ del Kalman è stimata dai residui di questo intervallo e resta poi fissata. Gli stimatori non ricevono i parametri veri della simulazione.

#### Lettura dei coefficienti

```{figure} ../assets/demo/demo_52_02.png
:alt: Coefficienti veri e stimati con finestre mobili, RLS e Kalman
:width: 100%
:name: fig-demo52-coefficienti

Demo 52: confronto dei coefficienti al ritardo 1. Il nero indica il valore vero. Il pannello $1\to3$ mostra le fluttuazioni attorno a una connessione diretta nulla.
```

Tutti i metodi ricostruiscono l'aumento e la diminuzione delle connessioni. Finestre e RLS oscillano intorno ai valori veri; con la taratura scelta, il Kalman produce traiettorie generalmente più regolari ma segue più lentamente parte delle transizioni. La variazione del coefficiente diagonale di $x_1$ riflette invece il cambiamento della sua frequenza propria.

Il pannello $1\to3$ mostra che un coefficiente stimato diverso da zero può comparire anche quando l'arco è assente nel modello generatore. La dimensione di queste fluttuazioni va valutata rispetto all'incertezza della stima.

#### Lettura delle mappe PDC

Per ogni modello locale si calcola la **PDC quadratica standard**:

$$
\operatorname{PDC}_{i\leftarrow j}(t,f)
=
\frac{|A_{ij}(f,t)|^2}
{\displaystyle\sum_{k=1}^{M}|A_{kj}(f,t)|^2}.
$$

La normalizzazione è per colonna, cioè per sorgente $j$, e comprende il termine diagonale. Il riferimento teorico si ottiene congelando i coefficienti veri al tempo considerato; le stime vengono valutate sulla stessa griglia tempo-frequenza.

```{figure} ../assets/demo/demo_52_03.png
:alt: PDC teorica e stimata per le direzioni 1 verso 2, 2 verso 3 e 1 verso 3
:width: 100%
:name: fig-demo52-pdc

Ogni riga rappresenta una direzione; le colonne mostrano il riferimento teorico e i tre stimatori. La scala cromatica è comune, da 0 a 1. Il tratteggio bianco indica la frequenza propria del primo oscillatore.
```

Nella prima riga l'interazione $1\to2$ si rafforza e successivamente si riduce; la frequenza variabile di $x_1$ modifica la regione in cui la PDC è elevata. Nella seconda riga diminuisce l'interazione $2\to3$. La terza rimane vicina a zero, coerentemente con l'assenza dell'arco diretto $1\to3$.

La larghezza delle regioni colorate dipende dalla dinamica autoregressiva e dalla normalizzazione. Il tratteggio indica la frequenza propria del modello generatore: non impone che ogni massimo della PDC debba coincidere esattamente con quella curva.

#### Verifica rispetto alla simulazione

Lo script controlla stabilità locale, normalizzazione della PDC, recupero dei livelli nei plateau ed errore rispetto ai parametri veri. Le metriche sono calcolate tra 15 e 115 s, dopo l'assestamento iniziale. Un unico VAR costante stimato sull'intera registrazione fornisce una baseline descrittiva; non è un metodo online.

| Modello | RMSE coefficienti | RMSE PDC fuori diagonale | Modelli locali stabili |
|---|---:|---:|---:|
| VAR costante | 0.0711 | 0.1706 | — |
| Finestre mobili | 0.0235 | 0.0402 | 100% |
| RLS | 0.0240 | 0.0371 | 100% |
| Kalman | 0.0209 | 0.0339 | 100% |

I valori si riferiscono alla realizzazione riproducibile verificata in MATLAB R2025b. I tre metodi adattivi descrivono questa dinamica meglio della baseline costante. Il risultato non stabilisce una superiorità generale del Kalman: rapidità delle transizioni, quantità di dati e taratura possono modificare le prestazioni relative.

```{admonition} Accuratezza della ricostruzione e significatività
:class: important
Le verifiche della Demo 52 confrontano le stime con un modello generatore noto. Non sono test di significatività. Le piccole PDC sulle connessioni assenti non devono essere interpretate come archi significativi senza una procedura inferenziale dedicata.
```

La [pagina della demo](../demo/demo_52_tv_var.md) contiene tutte le figure, i parametri e lo script scaricabile. Modificare un parametro alla volta fra $W$, $\lambda$ e $q$ permette di osservare il compromesso tra regolarità e adattamento.

## 17. Inferenza sulle mappe tempo-frequenza

Una mappa tempo-frequenza contiene un valore per molte combinazioni di direzione, tempo e frequenza. Prima di costruire una soglia bisogna definire la domanda statistica, perché ipotesi diverse richiedono distribuzioni nulle diverse.

| Domanda | Esempio di ipotesi nulla |
|---|---|
| esiste una connessione? | $\theta_{i\leftarrow j}(t,f)=0$ |
| la connessione cambia rispetto al basale? | $\theta_{i\leftarrow j}(t,f)-\theta_{i\leftarrow j}^{\mathrm{base}}(f)=0$ |
| due condizioni differiscono? | $\theta_{i\leftarrow j}^{(1)}(t,f)-\theta_{i\leftarrow j}^{(2)}(t,f)=0$ |

Un test di presenza della connessione non equivale a un test di variazione rispetto al basale. Una connessione può essere significativa in entrambe le condizioni senza che la loro differenza sia significativa; al contrario, due valori singolarmente incerti possono mostrare una differenza sistematica tra condizioni.

### 17.1 Surrogati e intera pipeline di stima

Quando si usano dati surrogati, ciascun dataset deve attraversare la stessa procedura applicata ai dati osservati:

1. preprocessing previsto dall'analisi;
2. stima del TV-VAR con gli stessi parametri adattivi;
3. controllo e gestione della stabilità;
4. trasformazione nel dominio della frequenza;
5. calcolo della misura e della statistica finale.

Applicare i surrogati soltanto alle mappe già calcolate ignora l'incertezza introdotta dalla stima dei coefficienti. Anche eventuali scelte guidate dai dati devono essere ripetute nei surrogati oppure fissate prima del test.

Il tipo di surrogato deve preservare gli aspetti che non appartengono all'ipotesi nulla. Per esempio, un test tra condizioni può usare permutazioni delle etichette quando il disegno sperimentale rende le etichette scambiabili. Un test di connessione richiede invece una procedura che elimini la relazione in esame conservando, per quanto possibile, la dinamica univariata e le altre caratteristiche rilevanti. La randomizzazione delle fasi descritta nella sezione 13 produce una null globale e può eliminare anche la struttura temporale legata all'evento; non è automaticamente adatta a ogni domanda tempo-frequenza.

### 17.2 Confronti multipli

Una soglia puntuale al livello $\alpha$ controlla l'errore di un singolo punto, ma non la probabilità di falsi positivi sull'intera mappa. Inoltre, punti vicini nel tempo e nella frequenza sono correlati per effetto della dinamica dei segnali e della memoria dello stimatore.

Con i surrogati si può usare una statistica del massimo. Per ogni replica $b$ si calcola

$$
T_b^{\max}
=
\max_{i\neq j,\,t,\,f}
\theta_{i\leftarrow j}^{*(b)}(t,f).
$$

Il quantile $1-\alpha$ della distribuzione di $T_b^{\max}$ fornisce una soglia che controlla l'errore familiare sulla famiglia di direzioni, tempi e frequenze inclusa nel massimo. La famiglia deve essere dichiarata: limitare il massimo a una banda o a una finestra temporale è corretto quando tale scelta è motivata a priori, non quando deriva dall'ispezione della mappa osservata.

I test basati su cluster sfruttano la contiguità dei punti. In ogni surrogato si applica una soglia iniziale, si identificano i cluster contigui e si registra la massima statistica di cluster, per esempio la somma dei valori sopra soglia. Il confronto riguarda quindi l'estensione o l'intensità complessiva della regione, non ciascun pixel separatamente. Il risultato consente di affermare che un cluster è significativo, ma non che ogni suo punto lo sia individualmente.

La FDR rappresenta un'altra possibilità quando l'obiettivo è controllare la proporzione attesa di falsi positivi tra i risultati dichiarati significativi. La scelta tra FDR, max-statistic e cluster dipende dalla domanda e dal tipo di conclusione che si vuole sostenere.

### 17.3 Interpretazione e analisi di sensibilità

È prudente evitare interpretazioni basate su singoli pixel isolati o su variazioni che coincidono con instabilità del modello, artefatti o grandi innovazioni. Una struttura fisiologicamente plausibile dovrebbe mostrare una certa coerenza temporale o spettrale e rimanere riconoscibile quando si modificano moderatamente i parametri dell'analisi.

L'analisi di sensibilità dovrebbe verificare almeno l'effetto di:

- ordine $p$ del VAR;
- parametro di dimenticanza o covarianza $\mathbf Q_{\theta}$;
- inizializzazione e intervallo di assestamento;
- preprocessing e selezione dei trial;
- famiglia di confronti e metodo di correzione.

Le bande di frequenza e le finestre temporali motivate dall'ipotesi sperimentale possono essere definite prima dell'analisi e usate per ottenere statistiche riassuntive. Questo riduce il numero di confronti e rende la conclusione più vicina alla domanda fisiologica, purché la definizione non venga adattata dopo avere osservato i risultati.

## 18. Checklist per un'analisi riproducibile

Un'analisi di connettività basata su VAR dovrebbe permettere a un altro ricercatore di ricostruire sia il modello sia tutte le decisioni che trasformano i segnali nelle mappe finali.

| Fase | Informazioni da riportare |
|---|---|
| Dati | segnali o sorgenti analizzati, frequenza di campionamento, durata, numero di soggetti e trial |
| Preprocessing | riferimento, filtri, rimozione di trend e artefatti, criteri di esclusione, eventuale standardizzazione |
| Struttura del VAR | variabili incluse, ordine $p$, intercetta, procedura di selezione e numero effettivo di osservazioni |
| Stima stazionaria | metodo di stima, regolarizzazione e trattamento dei segmenti o trial |
| Stima tempo-variante | finestre $W$ e passo $Q$, oppure $\lambda$, oppure modello di stato, inizializzazione, $\mathbf Q_{\theta}$ e $\mathbf R_u$ |
| Validazione | stabilità, test di Portmanteau sui residui, errori di previsione e intervalli esclusi perché inaffidabili |
| Connettività | definizione esatta di GC, PDC o DTF, convenzione degli indici, normalizzazione e griglia di frequenza |
| Inferenza | ipotesi nulla, procedura di generazione dei surrogati, numero di repliche, statistica e correzione multipla |
| Sensibilità | variazioni di ordine, parametri adattivi, preprocessing, banda e finestra temporale |
| Riproducibilità | software e versione, funzioni o codice utilizzati, seme casuale e disponibilità dei parametri completi |

Per i modelli tempo-varianti è inoltre opportuno dichiarare:

1. se le stime sono filtrate o smussate;
2. quale scala temporale effettiva deriva dai parametri adattivi;
3. come sono stati combinati eventuali trial;
4. come sono state gestite le stime localmente instabili;
5. quali dimensioni sono incluse nella famiglia dei confronti multipli.

```{admonition} Interpretazione finale
:class: warning
La causalità di Granger descrive informazione predittiva direzionale condizionata sulle variabili osservate e sul modello scelto. Una variazione della misura nel tempo non dimostra da sola una modifica di un meccanismo fisiologico. La conclusione richiede un modello locale adeguato, robustezza rispetto alle scelte analitiche e coerenza con il disegno sperimentale.
```
