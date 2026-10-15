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

## Sfruttare il seguento scheletro per completare l'esercitazione


```matlab
clear
clc

% simulazione u(t)
N = 100;
mu=[0 0];
sigma=[0.3 0; 0 0.3];
%potete provare a inserire un valore di cxy
%con |cxy|<sigmax*sigmay
u=mvnrnd(mu, sigma, N)';

x1(1)=0;%le condizioni iniziali possono essere diverse
x2(1)=0;%se avessi un modello di ordine 2 servirebbero altre due condizioni
X = [x1;x2];
A = [-0.7 -0.2; -0.2 -0.7];
p = 1;
for k=2:N
    X(:,k)=A*X(:,k-1)+u(:,k);
end

%da ora in poi consideriamo le due serie temporali e verifichiamo i valori degli indici di Granger
%dovremmo stimare due modelli marginali e uno completo


% =========================================================================
% 1. MODELLO COMPLETO (VAR)
% =========================================================================
Y = X(:, p+1:end);
[M, T] = size(X);
num_obs = T - p;

% Costruzione di Z per il modello completo
Z = zeros(M*p + 1, num_obs);
Z(1, :) = 1;
for k = 1:p
    row_idx = (1 + (k-1)*M + 1) : (1 + k*M);
    Z(row_idx, :) = X(:, (p - k + 1):(T - k));
end

B = Y * Z' / (Z * Z'); % Nota: / è preferibile numericamente rispetto a *inv()
u_hat = Y - B * Z;

% Varianze dei residui del modello completo (sulla diagonale della covariance)
sigma1_full = var(u_hat(1, :));
sigma2_full = var(u_hat(2, :));


% =========================================================================
% 2. MODELLI MARGINALI (AR UNIVARIATI)
% =========================================================================

```

## Domande per la relazione

1. Quale struttura direzionale ci si aspetta dai tre sistemi?
2. Quanto cambia l'errore $\|\widehat{\mathbf{A}}_1-\mathbf{A}_1\|_F$ passando da 100 a 1000 campioni?
3. Un valore di GC piccolo può essere significativo con $N$ grande? Spiegare.
4. Il mancato rifiuto dell'ipotesi nulla dimostra l'assenza di connessione?
5. Quale modello nullo viene usato dal test scelto?

## Extra
Prendere spunto  dalla [Demo 45](../demo/demo_45_validazione.md) per validare il modello. 

