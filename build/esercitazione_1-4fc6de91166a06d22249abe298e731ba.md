# Esercitazione 1

## Simulazione e stima di un VAR(1) bivariato

### Obiettivi

L'esercitazione consolida simulazione, stima ai minimi quadrati, controllo della stabilità e causalità di Granger nel dominio del tempo. Il confronto tra tre matrici permette di collegare la struttura dei coefficienti alla direzione attesa della GC.

## Modelli

Per ciascuno dei tre sistemi usare

$$
\mathbf{x}_t=\mathbf{A}_1\mathbf{x}_{t-1}+\mathbf{u}_t,
\qquad
\mathbf{u}_t\sim\mathcal{N}(\mathbf{0},\boldsymbol{\Sigma}),
$$

con

$$
\boldsymbol{\Sigma}=
\begin{bmatrix}
0.3&0\\
0&0.3
\end{bmatrix}.
$$

### Sistema A

$$
\mathbf{A}_1=
\begin{bmatrix}
-0.7&0\\
-0.2&-0.7
\end{bmatrix}.
$$

### Sistema B

$$
\mathbf{A}_1=
\begin{bmatrix}
-0.7&-0.2\\
-0.2&-0.7
\end{bmatrix}.
$$

### Sistema C

$$
\mathbf{A}_1=
\begin{bmatrix}
-0.7&-0.2\\
0&-0.7
\end{bmatrix}.
$$

Ricordare che $A_{ij}$ descrive $j\rightarrow i$. Prima di simulare, formulare una previsione qualitativa sulla GC nei tre sistemi.

## Attività richieste

Per ogni sistema:

1. verificare la stabilità mediante gli autovalori di $\mathbf{A}_1$;
2. simulare le innovazioni gaussiane con covarianza $\boldsymbol{\Sigma}$;
3. generare il processo ricorsivo;
4. rappresentare $x_1(t)$ e $x_2(t)$;
5. stimare un VAR(1) dai dati simulati;
6. confrontare $\widehat{\mathbf{A}}_1$ con la matrice originale;
7. stimare la covarianza residua $\widehat{\boldsymbol{\Sigma}}$;
8. calcolare la GC nelle due direzioni;
9. valutarne la significatività statistica;
10. confrontare i risultati dei tre sistemi.

Usare $\Delta t=10$ ms e ripetere con:

- $N=100$ campioni;
- $N=1000$ campioni.

Valutare l'effetto di $N$ su $\widehat{\mathbf{A}}_1$, variabilità delle stime e direzione della GC. Generare 200-500 campioni iniziali da scartare come burn-in.

## Modalità di lavoro

Scrivere autonomamente uno script MATLAB che svolga nell'ordine le attività richieste. Il codice di partenza e la soluzione non sono forniti in questa versione della dispensa. Motivare nel codice o nella relazione le scelte relative a inizializzazione, burn-in, costruzione dei regressori, stima dei modelli ridotti e test di significatività.

## Codice MATLAB da condividere

In questo spazio verrà inserito il solo frammento di codice MATLAB che sarà reso disponibile agli studenti.

```matlab
% Codice MATLAB da inserire
```

## Domande per la relazione

1. Quale struttura direzionale ci si aspetta dai tre sistemi?
2. Quanto cambia l'errore $\|\widehat{\mathbf{A}}_1-\mathbf{A}_1\|_F$ passando da 100 a 1000 campioni?
3. Un valore di GC piccolo può essere significativo con $N$ grande? Spiegare.
4. Il mancato rifiuto dell'ipotesi nulla dimostra l'assenza di connessione?
5. Quale modello nullo viene usato dal test scelto?

## Consegna minima

La relazione deve contenere matrici vere e stimate, autovalori, covarianza residua, GC nelle due direzioni, esito dell'inferenza e un breve commento sull'effetto della numerosità.

