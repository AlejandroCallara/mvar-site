# Esercitazione 2

## PDC e DTF in un sistema VAR multivariato

### Obiettivo

Simulare un VAR(2) a quattro variabili, stimarne congiuntamente i coefficienti e calcolare PDC, gPDC, DTF e dDTF nel dominio della frequenza.

## Modello

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

Scrivere autonomamente uno script MATLAB che svolga nell'ordine le attività richieste. Il codice di partenza e la soluzione non sono forniti in questa versione della dispensa. Motivare nel codice o nella relazione le scelte relative a simulazione, stima congiunta, griglia delle frequenze, normalizzazione delle misure e rappresentazione grafica.

## Domande per la relazione

1. Quali connessioni sono rappresentate direttamente da coefficienti fuori diagonale?
2. Perché una DTF $1\rightarrow3$ può essere diversa da zero anche senza coefficiente diretto?
3. Quale normalizzazione usa la PDC e quale usa la DTF?
4. Perché PDC e gPDC coincidono con innovazioni a varianza unitaria?
5. Quale effetto produce una covarianza delle innovazioni non diagonale?

## Consegna minima

Includere matrici vere e stimate, raggio spettrale, covarianza residua, grafici PDC/DTF/dDTF, confronto PDC-gPDC e interpretazione della rete.

