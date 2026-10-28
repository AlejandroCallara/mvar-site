# Metodi per l'analisi di segnali multidimensionali

## VAR e causalità di Granger

**Laurea Magistrale in Ingegneria Biomedica**  
**Università di Pisa**  
**Docenti:** Dr. Alejandro Luis Callara e Prof. Nicola Vanello  
**Durata:** 12 ore, con teoria, simulazioni ed esercitazioni MATLAB

Questa dispensa accompagna il modulo dedicato ai modelli vettoriali autoregressivi e alla causalità di Granger. Il percorso parte dal significato predittivo della causalità di Granger, introduce la stima e la validazione dei modelli VAR, passa alle misure nel dominio della frequenza e conclude con i modelli tempo-varianti.

Il corso usa MATLAB. Le demo sono fornite come script `.m` leggibili e scaricabili. Ogni pagina descrive lo scopo della simulazione, i risultati attesi e le cautele interpretative.

## Obiettivi formativi

Al termine del modulo lo studente saprà:

- formulare e stimare un modello VAR multivariato;
- scegliere l'ordine del modello e verificarne stabilità e adeguatezza;
- calcolare e interpretare la causalità di Granger nel dominio del tempo;
- distinguere dipendenza istantanea, influenza direzionale e causalità fisica;
- ricavare la rappresentazione spettrale del VAR e calcolare PDC e DTF;
- costruire procedure di inferenza con test parametrici, bootstrap e dati surrogati;
- riconoscere gli effetti di SNR, numerosità campionaria e variabili non osservate;
- confrontare finestre mobili, RLS e filtro di Kalman su un TV-VAR e interpretare la PDC tempo-frequenza.

## Percorso consigliato

1. Leggere il capitolo principale seguendo l'ordine delle sezioni.
2. Eseguire le demo MATLAB quando vengono richiamate nel testo.
3. Svolgere le tre esercitazioni: [GC nel tempo](esercitazioni/esercitazione_1.md), [misure nel dominio della frequenza](esercitazioni/esercitazione_2.md) e [inferenza mediante surrogati di fase](esercitazioni/esercitazione_3.md).

**Le esercitazioni saranno argomento di esame.**

```{admonition} Convenzione degli indici
:class: important
In questa dispensa, nelle matrici VAR e nelle misure spettrali, l'elemento `(i,j)` descrive la direzione dalla sorgente `j` al bersaglio `i`.
```

