# Esercitazione 3 — Significatività di PDC e DTF mediante surrogati di fase

## Obiettivo

Completare l'analisi dell'[esercitazione 2](esercitazione_2.md) costruendo una distribuzione nulla per PDC e DTF mediante randomizzazione indipendente delle fasi dei canali. Scrivere autonomamente uno script MATLAB che generi i surrogati, ristimi il modello e confronti le misure osservate con soglie empiriche.

Mentre la seconda esercitazione riguardava simulazione, stima e interpretazione delle misure, questa riguarda l'inferenza statistica. 

Per le misure nel dominio della frequenza e la loro interpretazione:

Per l'ipotesi nulla dei surrogati, le soglie e i confronti multipli, fare riferimento alla sezione 13 della dispensa.

## Dati e impostazioni

Riprendere il VAR(2) a quattro variabili dell'esercitazione 2, nel caso iniziale con innovazioni gaussiane indipendenti e covarianza $\boldsymbol\Sigma_u=\mathbf I_4$. Usare una singola realizzazione organizzata come matrice $M\times N$, con $M=4$, $N=1000$ campioni dopo il transitorio e $\Delta t=1$ s.

Riutilizzare la stessa realizzazione e le stime ottenute nell'esercitazione 2. Usare $p=2$ e una griglia di 128 frequenze tra zero e Nyquist, mantenuta identica per dati originali e surrogati.

Le misure richieste sono le versioni quadratiche di PDC e DTF.

## Parte A — Misure osservate, ipotesi nulla e generazione dei surrogati

1. Calcolare PDC e DTF dal VAR stimato sui dati originali. Conservare le misure osservate come riferimento per tutti i confronti successivi.
2. Esplicitare l'ipotesi nulla della randomizzazione indipendente delle fasi: si conservano le proprietà spettrali univariate, ma si distruggono le relazioni sistematiche di fase fra canali.
3. Spiegare perché questa è una null globale di assenza di dipendenza tra canali e non una null selettiva che elimina soltanto un arco mantenendo inalterato il resto della rete.
4. Per ogni canale, calcolare la FFT e costruire nuovi coefficienti con lo stesso modulo e fasi casuali indipendenti. Randomizzare soltanto le fasi delle componenti a frequenza positiva strettamente comprese tra zero e $(f_s/2)$. Lasciare invariati i coefficienti a frequenza zero e a $(f_s/2)$, e costruire le componenti a frequenza negativa mediante simmetria coniugata. Se si sottrae la media prima della FFT, riaggiungerla dopo la trasformata inversa.
5. Verificare su almeno un surrogato che il modulo FFT coincida con quello originale. Confrontare segnali e periodogrammi; verificare anche media e varianza.
6. Generare inizialmente $B=200$ dataset surrogati per controllare lo script, poi aumentare a $B=1000$ per l'analisi finale. Le fasi devono essere indipendenti fra canali e fra repliche; non reimpostare lo stesso seme randomico all'interno del ciclo.

## Parte B — Ristima e distribuzioni nulle

7. Per ogni dataset surrogato, ristimare il VAR(2) con la stessa procedura usata sui dati originali.
8. Ricalcolare entrambe le misure sulla stessa griglia di frequenze e raccogliere i risultati per direzione, frequenza e replica. 
9. Per la direzione $1\to2$, usare il bin della griglia più vicino a 0.1 Hz: con 128 frequenze uniformi tra 0 e 0.5 Hz è il bin a circa 0.0984 Hz. Rappresentare due istogrammi: uno dei valori di PDC ottenuti dai surrogati e uno dei valori di DTF. Su ciascun istogramma indicare, con una linea verticale, il corrispondente valore calcolato sui dati originali. Commentare la posizione del valore osservato rispetto alla distribuzione surrogata.

## Parte C — Soglie puntuali e risultati

10. Fissare $(\alpha=0.05)$. Per ogni direzione fuori diagonale e ogni frequenza, ordinare in senso crescente i $(B)$ valori ottenuti dai surrogati e usare come soglia il valore in posizione $(k=\lceil0.95B\rceil)$. Calcolare separatamente le soglie per PDC e DTF e confrontarle con i rispettivi valori osservati. Notare che le soglie sono puntuali e non correggono i confronti multipli fra direzioni e frequenze.
11. Calcolare anche il p-value empirico monolaterale:

    $$
    p_{i\leftarrow j}(f)=\frac{1+\sum_{b=1}^{B}\mathbb I[\theta^{*(b)}_{i\leftarrow j}(f)\geq\theta^{\mathrm{obs}}_{i\leftarrow j}(f)]}{B+1},
    $$
    
    dove $\theta$ indica, separatamente, PDC quadratica o DTF quadratica. 
12. Produrre due matrici $4\times4$ di grafici: PDC osservata con soglia e DTF osservata con soglia. 
13. Discutere i risultati rispetto alla rete assegnata, distinguendo connessioni dirette e propagazioni indirette. Ricordare che una DTF sopra soglia non dimostra l'esistenza di un arco diretto; un mancato superamento non dimostra l'assenza di una connessione.

```{admonition} Interpretazione delle soglie
:class: important
Le soglie richieste sono puntuali: riguardano una misura, una direzione e una frequenza alla volta. Non controllano al 5% l'errore sull'intera matrice e su tutta la griglia. Inoltre i surrogati distruggono tutte le relazioni fra canali: il confronto non è automaticamente un test condizionale del singolo arco in una rete con altre connessioni presenti.
```

## Soluzione

### Impostazioni dell'esecuzione di riferimento

I risultati seguenti sono basati sul modello dell'esercitazione 2, con $N=1000$ campioni dopo un transitorio di 500 campioni, $p=2$, $\Delta t=1$ s, innovazioni indipendenti a varianza unitaria e seme casuale 11. Per ciascuno dei $B=1000$ surrogati si ristima il VAR e si ricalcolano PDC e DTF quadratiche sulle stesse 128 frequenze uniformemente distribuite tra 0 e 0.5 Hz.

Il raggio spettrale della matrice compagna è 0.9065 per il modello vero e 0.9168 per quello stimato: entrambi risultano stabili. I valori riportati illustrano una singola realizzazione, non risultati numerici che debbano essere identici per qualunque simulazione. Per un confronto esatto occorre usare gli stessi dati, lo stesso generatore casuale e la stessa sequenza di operazioni.

### Verifica delle proprietà univariate

La randomizzazione cambia la forma d'onda dei segnali, ma mantiene i moduli dei coefficienti FFT di ciascun canale. Di conseguenza si conservano il periodogramma, la media e la varianza. Le relazioni sistematiche fra le fasi dei diversi canali vengono invece distrutte.

![Confronto fra un segnale originale e un surrogato e sovrapposizione dei periodogrammi](../assets/esercitazione_3/verifica_surrogati.png)

La figura mostra il primo canale e l'ultimo surrogato generato. Le due curve spettrali coincidono entro la precisione numerica; la sovrapposizione non indica che i segnali nel tempo siano uguali.

| Verifica sui quattro canali dell'ultimo surrogato | Errore massimo |
|---|---:|
| Moduli FFT, errore assoluto massimo diviso per il massimo modulo originale | $3.83\cdot10^{-16}$ |
| Media, differenza assoluta | $1.18\cdot10^{-16}$ |
| Varianza, differenza relativa per canale | $2.30\cdot10^{-15}$ |


### Distribuzioni nulle: direzione $1\to2$ a 0.0984 Hz

![Distribuzioni surrogate della PDC e della DTF per la direzione 1 verso 2 al bin di 0.0984 Hz](../assets/esercitazione_3/distribuzioni_nulle.png)

| Misura, direzione $1\to2$, frequenza 0.0984 Hz | Valore osservato | Soglia, posizione 950 | P-value empirico |
|---|---:|---:|---:|
| PDC quadratica | 0.5215 | 0.4571 | $27/1001\approx0.0270$ |
| DTF quadratica | 0.9874 | 0.5488 | $1/1001\approx0.0010$ |


### PDC e DTF osservate e soglie puntuali

Nelle matrici seguenti, il pannello $(i,j)$ rappresenta $j\to i$. La curva blu è la misura osservata; il tratteggio nero è la soglia empirica in posizione $\lceil0.95B\rceil=950$. La diagonale non viene sottoposta al test sui collegamenti, ma resta inclusa nei denominatori delle normalizzazioni.

![PDC quadratica osservata e soglie puntuali nelle dodici direzioni fuori diagonale](../assets/esercitazione_3/risultati_01.png)

![DTF quadratica osservata e soglie puntuali nelle dodici direzioni fuori diagonale](../assets/esercitazione_3/risultati_02.png)

La PDC individua contributi diretti coerenti con i coefficienti del modello: $1\to2$, $3\to2$, $2\to3$ e $1\to4$. Il contributo $2\to3$ dipende dalla frequenza anche per la presenza di due ritardi, che possono attenuarsi reciprocamente in alcune bande.

La DTF comprende invece la propagazione attraverso tutti i cammini. Il caso più evidente è $1\to3$: il coefficiente diretto è assente, ma l'influenza del nodo 1 raggiunge il nodo 3 attraverso il nodo 2. Una DTF elevata o sopra soglia su questa direzione non richiede quindi di aggiungere un arco diretto al grafo.



## Approfondimento facoltativo

Per ogni surrogato, individuare il valore massimo della PDC considerando tutte le direzioni tra nodi diversi e tutte le frequenze. Usare il 95° percentile dei massimi ottenuti come un’unica soglia per l’intera matrice PDC. Ripetere separatamente la procedura per la DTF.
Confrontare queste soglie con quelle puntuali e osservare come cambiano i risultati. Le soglie puntuali riguardano una direzione e una frequenza alla volta; la soglia del massimo considera insieme tutte le direzioni e frequenze della stessa misura.
Attenzione: i surrogati rappresentano una rete priva di relazioni fra canali. La soglia così costruita controlla i superamenti rispetto a questa ipotesi globale; non garantisce automaticamente lo stesso controllo sui singoli collegamenti assenti in una rete che contiene altre connessioni.
