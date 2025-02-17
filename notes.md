# IP Notizen 

## Literaturrecherche 

### Alternierende Optimierungsprobleme

+ Schwierige Optimierung über alle $s$ Variablen 
+ Eingangsfunktion hat mehrerer Veränderliche

$$f(x_1, x_2, ..., x_n)$$

## Formulierung unseres Ansatzes:

### Voraussetzungen 

+ $m$ Normalverteile, voneinander unabhängige Zufallsvariablen $$X$$
mit Realisierungen $x_1, x_2, ..., x_n$
+ Anzahl der Summanden wird o.B.d.A. auf $n=3$ festgesetzt. 
+ Eine Linearkombination kann dann definniert werden als $$f_1 = h_1 \cdot x_1 + h_2 \cdot x_2 + h_3 \cdot x_3$$
$f_1$ definiert hier dann die erste aus $\lfloor \frac{m}{n}\rfloor$ Linearkombinationen
+ Definiere die Grenzen der Quantisierer als $j$ Grenzen im Vektor $\mathbf{b}$ mit $|\mathbf{b}| = j$
+ Abhängngig von der Anzahl der Summanden $n$ ergeben sich bei binären Weights $2^m$ mögliche Linearkombinationen
+ Die Funktionen $f_i$ bilden ein Skalarfeld mit einer Abbildung $$\real^n \rightarrow \real$$

### Definition Optimierungsproblem

$$argmax_h min_j|\mathbf{h}^T\mathbf{x} -b_j| \text{ s.t.} h \in {\pm 1}$$

-> Finde die Kombination an Weights, so dass der Abstand zur nächstgelegenen Grenze maximiert wird.

Grenzen $b_j$ können nicht von anfang an bekannt sein, müssen geraten werden über $n$-Fache Faltung der Eingangsverteilung mit sich selbst

Problem: Wahl der Weights nach dieser Vorschrift konvergiert nicht sinnvoll in eine neue Verteilung

-> Nebenbedingung gefordert

Mögliche Optionen: 

+ Speichere Weights für gleiche Abstände zu Grenzen ab und balance regelmäßig die Grenzen wieder aus je mehr man dazu nimmmt 

Ziel: Erreiche möglichst eine Gleichvertilung der Weights an die ausgangsverteilung über Normalverteilungen

## Aktuelle Punkte 

Erste Bach Iteration resultiert in einer sehr ungleichmäßigen Quantisierung

- Man könnte hier noch dran arbeiten und diese Verteilung gleichmäßiger basteln, ohne alternierend das zu verbessern 
- Iterativ versuchen die Grenzen neu zu sezten, damit es dann am Ende doch einfach passt. 
- 