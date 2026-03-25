## Méthodes : Reconstruction conjointe du champ de conduction et du cycle organisateur

### Cadre

Soit \(M\) une surface triangulée 2D immergée en 3D, représentant le myocarde, munie de sa métrique induite \(g_0\).  
On observe une carte d’activation \(\tau_{\mathrm{obs}} : M \to \mathbb{R} / T\mathbb{Z}\) (temps modulo la période \(T\)).

On suppose que la conduction est **isotrope mais hétérogène**, décrite par une vitesse scalaire
\[
c : M \to (0,+\infty),
\]
induisant une métrique électrophysiologique
\[
g_c = c^{-2} g_0.
\]

Dans ce cadre, les temps d’activation idéaux satisfont une équation eikonale sur la surface :
\[
|\nabla_M \tau(x)|_{g_0} = \frac{1}{c(x)}.
\]

Le cycle organisateur \(\gamma\) est modélisé comme une **géodésique fermée** de la métrique \(g_c\), de longueur
\[
L_{g_c}(\gamma) = \int_\gamma \frac{ds}{c(x)}.
\]

---

### Variables inconnues

On cherche à estimer conjointement :

- le champ de conduction \(c(x)\),
- une boucle fermée \(\gamma \subset M\), dans une classe d’homotopie donnée.

Pour des raisons numériques, on paramètre
\[
u(x) = \log c(x), \quad c(x)=e^{u(x)}.
\]

---

### Prétraitement : phase et déroulement du temps

On définit la phase complexe :
\[
\phi(x) = e^{i 2\pi \tau_{\mathrm{obs}}(x)/T}.
\]

La phase est utilisée pour :

- détecter les singularités (défauts de phase),
- identifier des régions candidates pour le cœur organisateur,
- sélectionner un ensemble de classes d’homotopie plausibles.

Un **déroulement de phase** (phase unwrapping) est ensuite effectué pour reconstruire un temps \(\tau(x)\) continu sur \(M\) (à une constante additive près), en évitant les zones singulières.

---

### Initialisation

#### Boucles initiales

Pour chaque classe topologique \(\alpha\) jugée pertinente, on calcule une géodésique fermée minimale \(\gamma_\alpha^{(0)}\) pour la métrique \(g_0\) (i.e. \(c \equiv 1\)), en utilisant des méthodes de type Crane et al.

Ces boucles servent d’initialisations multiples.

#### Champ de conduction initial

On estime un champ initial :
\[
c_0(x) \approx \frac{1}{|\nabla_M \tau(x)|_{g_0}},
\]
calculé par éléments finis (P1) sur le maillage.

On projette ensuite sur des bornes physiologiques :
\[
c_{\min} \le c_0(x) \le c_{\max},
\]
et on lisse ce champ (Laplacien de surface).

---

### Formulation variationnelle

On définit une fonctionnelle d’énergie jointe :
\[
\mathcal J(u,\gamma)
=
\mathcal J_{\mathrm{data}}(u,\gamma)
+
\lambda \mathcal J_{\mathrm{reg}}(u)
+
\mu \mathcal J_{\mathrm{geo}}(u,\gamma)
+
\nu \mathcal J_{\mathrm{phys}}(u,\gamma).
\]

#### 1. Terme d’adéquation aux données

On impose la cohérence avec l’équation eikonale :
\[
\mathcal J_{\mathrm{data}}(u)
=
\int_M
\left(
|\nabla_M \tau(x)|_{g_0} - e^{-u(x)}
\right)^2 dA.
\]

Ce terme est évalué uniquement sur les régions où \(\tau\) est régulière (exclusion des singularités et collisions de fronts).

---

#### 2. Régularisation spatiale

On pénalise les variations de \(u\) :
\[
\mathcal J_{\mathrm{reg}}(u)
=
\int_M |\nabla_M u|^2 dA.
\]

Sur le mesh, cela correspond à :
\[
u^\top L u,
\]
où \(L\) est le Laplacien cotangent.

---

#### 3. Terme géométrique (cycle organisateur)

On impose que \(\gamma\) soit une géodésique fermée de \(g_c\) et qu’elle respecte la période :

\[
\mathcal J_{\mathrm{geo}}(u,\gamma)
=
\underbrace{\int_\gamma e^{-u(x)} ds}_{L_{g_c}(\gamma)}
+
\beta \, \mathcal J_{\mathrm{closure}}(\gamma)
+
\eta \, \mathcal J_{\mathrm{geod}}(u,\gamma).
\]

- \(\mathcal J_{\mathrm{closure}}\) pénalise la non-fermeture,
- \(\mathcal J_{\mathrm{geod}}\) pénalise l’écart à l’équation géodésique (courbure géodésique non nulle),
- on peut imposer explicitement :
\[
L_{g_c}(\gamma) \approx T.
\]

---

#### 4. Contraintes physiologiques

On impose :

- bornes :
\[
u_{\min} \le u(x) \le u_{\max},
\]
- éventuellement une pénalisation des gradients trop forts,
- exclusion des zones non conductrices (scar).

---

### Schéma d’optimisation

On utilise une **optimisation alternée** :

1. **Initialisation** : \(u^{(0)}, \gamma^{(0)}_\alpha\)

2. Pour chaque classe \(\alpha\), itérer :

   **(a) Mise à jour du champ \(u\)**  
   Minimiser \(\mathcal J(u,\gamma^{(k)})\) avec \(\gamma\) fixé  
   (descente de gradient, L-BFGS, ou Gauss–Newton).

   **(b) Mise à jour de la boucle \(\gamma\)**  
   Calculer une géodésique fermée dans la métrique
   \[
   g_c = e^{-2u^{(k+1)}} g_0
   \]
   dans la classe d’homotopie \(\alpha\).

3. Critère d’arrêt :
\[
\|u^{(k+1)} - u^{(k)}\| < \varepsilon,
\quad
d(\gamma^{(k+1)},\gamma^{(k)}) < \varepsilon.
\]

---

### Sélection de solution

On obtient une famille de solutions \((u_\alpha,\gamma_\alpha)\).

On sélectionne la solution finale selon :

- l’énergie totale minimale,
- la cohérence avec la période \(T\),
- la plausibilité physiologique (valeurs de \(c\)),
- la stabilité vis-à-vis des perturbations des données.

---

### Remarques

- Le problème est **non unique** à partir d’une seule carte d’activation : la régularisation est essentielle.
- La phase est utilisée pour la **topologie et l’initialisation**, mais l’inversion est effectuée sur le temps déroulé.
- La paramétrisation \(u=\log c\) améliore la stabilité numérique et garantit \(c>0\).