# Esercitazione 2 - Stima di misure di connettività nel dominio della frequenza

## PDC e DTF in un sistema VAR multivariato

## Obiettivo

Simulare un VAR(2) a quattro variabili, stimarne congiuntamente i coefficienti e calcolare PDC, gPDC, DTF e dDTF nel dominio della frequenza.

Questo stesso modello verrà utilizzato per affrontare l'analisi di significatività mediante surrogati nell'[esercitazione 3](esercitazione_3.md).


## Modello: esempio VAR(2) multivariato

Il sistema a quattro processi riprende l'esempio multivariato di Faes e Nollo (2011).

Generare $N=1000$ campioni con $\Delta t=1$ s:

$$
\mathbf{x}(t)=\mathbf{A}_1\mathbf{x}(t-1)+
\mathbf{A}_2\mathbf{x}(t-2)+\mathbf{u}(t),
$$

$$
\mathbf{u}(t)\sim\mathcal{N}(\mathbf{0},\mathbf{I}).
$$

Le equazioni sono

$$
\begin{aligned}
x_1(t)&=2\rho\cos(2\pi f_1)x_1(t-1)-\rho^2x_1(t-2)+u_1(t),\\
x_2(t)&=x_1(t-1)+0.5x_3(t-1)+u_2(t),\\
x_3(t)&=2\rho\cos(2\pi f_3)x_3(t-1)-\rho^2x_3(t-2)\\
&\quad+0.5x_2(t-1)+0.5x_2(t-2)+u_3(t),\\
x_4(t)&=x_1(t-2)+u_4(t),
\end{aligned}
$$

con

$$
\rho=0.9,\qquad f_1=0.1\ \mathrm{Hz},\qquad f_3=0.3\ \mathrm{Hz}.
$$

## Attività richieste

1. Generare le quattro serie eliminando un transitorio iniziale.
2. Costruire

   $$
   \mathbf{Y}=
   \begin{bmatrix}
   \mathbf{x}(p+1)&\cdots&\mathbf{x}(N)
   \end{bmatrix},
   $$

   $$
   \mathbf{Z}=
   \begin{bmatrix}
   \mathbf{x}(p)&\cdots&\mathbf{x}(N-1)\\
   \mathbf{x}(p-1)&\cdots&\mathbf{x}(N-2)
   \end{bmatrix}.
   $$

3. Stimare

   $$
   \widehat{\mathbf{B}}=
   \begin{bmatrix}
   \widehat{\mathbf{A}}_1&\widehat{\mathbf{A}}_2
   \end{bmatrix}
   =\mathbf{Y}\mathbf{Z}^{\mathsf T}
   (\mathbf{Z}\mathbf{Z}^{\mathsf T})^{-1}.
   $$

4. Confrontare le matrici stimate con quelle vere.
5. Costruire

   $$
   \mathbf{A}(f)=\mathbf{I}-\sum_{r=1}^{2}
   \widehat{\mathbf{A}}_re^{-i2\pi fr\Delta t},
   \qquad
   \mathbf{H}(f)=\mathbf{A}(f)^{-1}.
   $$

6. Calcolare

   $$
   PDC_{i\leftarrow j}^2(f)=
   \frac{|A_{ij}(f)|^2}{\sum_{k=1}^{M}|A_{kj}(f)|^2},
   $$

   $$
   DTF_{i\leftarrow j}^2(f)=
   \frac{|H_{ij}(f)|^2}{\sum_{k=1}^{M}|H_{ik}(f)|^2}.
   $$

7. Rappresentare PDC e DTF con due matrici $4\times4$ di grafici. Il pannello $(i,j)$ indica $j\rightarrow i$.
8. Discutere connessioni dirette e propagazioni indirette.
9. Calcolare la dDTF e verificare se l'eventuale trasferimento indiretto $1\rightarrow3$ viene attenuato.
10. Ripetere con

    $$
    \boldsymbol{\Sigma}_u=\operatorname{diag}(1,4,0.25,2)
    $$

    e confrontare PDC e gPDC.

## Modalità di lavoro

Scrivere autonomamente uno script MATLAB che svolga nell'ordine le attività richieste. Motivare nel codice o nella relazione le scelte relative a simulazione, stima congiunta, griglia delle frequenze, normalizzazione delle misure e rappresentazione grafica.


## Alcune considerazioni

La struttura contiene $1\rightarrow2$, $3\rightarrow2$, $2\rightarrow3$ e $1\rightarrow4$. Il percorso $1\rightarrow2\rightarrow3$ può emergere nella DTF anche quando non esiste un coefficiente diretto $1\rightarrow3$.

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

## Prosecuzione: inferenza statistica

L'[esercitazione 3](esercitazione_3.md) riprende lo stesso modello e la stessa realizzazione per costruire surrogati a fase randomizzata, ristimare il VAR e ricavare soglie puntuali per PDC e DTF. Le figure di confronto di questa pagina illustrano le misure, non sostituiscono l'analisi della loro significatività.

## Riferimento

Faes, L., & Nollo, G. (2011). Multivariate frequency domain analysis of causal interactions in physiological time series. In *Biomedical Engineering, Trends in Electronics, Communications and Software*.


