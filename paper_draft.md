# Working version

# Organizing Cycles of Cardiac Reentry as Closed Geodesics of an Electrophysiological Metric

## 1. Introduction

Cardiac reentry is a fundamental mechanism underlying many arrhythmias. In anatomically anchored reentry, electrical activation propagates periodically along a circuit constrained by the geometry and electrophysiological properties of the tissue.

A stable reentrant wave does not correspond to a single closed wavefront. Instead, it generates a **periodic ribbon of wavefronts** surrounding a central geometric structure. 

The central hypothesis of this work is that the **organizing cycle** of such a ribbon can be modeled as a **closed geodesic of an electrophysiological metric** defined on the myocardial domain.

This geometric viewpoint leads naturally to a variational description of reentrant circuits and connects the structure of possible circuits to the topology of the tissue domain.

## 2. Electrophysiological Metric

In the eikonal approximation of cardiac excitation, the propagation of activation fronts can be described in terms of local conduction velocities. In anisotropic cardiac tissue, conduction velocity depends on both spatial location and propagation direction, reflecting the underlying fiber architecture and electrophysiological heterogeneities.

Following the geometric formulation proposed by Young and Panfilov, anisotropic conduction can be naturally interpreted in terms of a **Riemannian metric** on the excitable domain.

Let

$$
(M,g)
$$

be a compact Riemannian manifold representing the myocardial domain.

* $M$ denotes the excitable tissue region.
* $g$ is the **electrophysiological metric**, encoding the anisotropic conduction properties of the tissue.

### 2.1 Conduction Velocity and Metric Tensor

In anisotropic tissue, the local conduction velocity can be described by a symmetric positive-definite tensor field

$$
D(x),
$$

often referred to as the **diffusion tensor** in reaction–diffusion models of cardiac excitation.

In the eikonal limit, wavefront propagation is governed by an anisotropic eikonal equation of the form

$$
\nabla \tau(x)^T D(x) \nabla \tau(x) = 1,
$$

where $\tau(x)$ denotes the activation time.

Young and Panfilov observed that this equation can be interpreted geometrically by introducing the metric tensor

$$
g = D^{-1}.
$$

thus the metric is: 

$$
g(\tau,\tau)=\tau^TD^{−1}\tau
$$



With this identification, the eikonal equation becomes

$$
g^{-1}(d\tau,d\tau) = 1.
$$

Thus activation time behaves as a **distance function in the metric&#x20;****$g$**.

This observation provides a natural geometric framework for describing anisotropic cardiac conduction.

### 2.2 Electrophysiological Length

For a smooth closed curve

$$
\gamma : S^1 \to M
$$

we define its electrophysiological length

$$
L_g(\gamma) =
\int_0^1 \|\dot{\gamma}(t)\|_g \, dt.
$$

Physiologically, this quantity represents the **conduction time required for an activation front to traverse the&#x20;**~~**path&#x20;**~~**&#x20;curve****$\gamma$**.

Minimizing curves for this functional correspond to **geodesics of the electrophysiological metric**, which represent locally optimal propagation trajectories.

### 2.3 Relation to Wavefront Propagation

In the eikonal regime, activation fronts propagate orthogonally to the gradient of activation time and follow paths that locally minimize travel time.

Consequently, the electrophysiological metric provides a natural geometric structure for describing the spatial organization of propagation.

In particular:

* shortest paths in the metric $g$ correspond to fastest conduction trajectories,
* anisotropy of conduction is encoded directly in the metric tensor.

This geometric formulation provides the basis for interpreting organizing cycles of reentrant waves as **closed geodesics of the electrophysiological metric**.

## 3. Free Loop Space

Let

$$
\Lambda M = H^1(S^1,M)
$$

denote the **free loop space** of $M$, consisting of $H^1$ Sobolev loops.

The space decomposes into connected components indexed by free homotopy classes:

$$
\Lambda M = \bigsqcup_{\alpha \in [S^1,M]} \Lambda_\alpha M .
$$

Each component $\Lambda_\alpha M$ corresponds to a distinct **topological class of circuits**.

This decomposition plays a central role in classifying potential organizing cycles of reentrant waves.

## 4. Periodic Reentrant Waves

A stable reentry produces a family of wavefronts

$$
\Gamma_t
$$

with temporal period

$$
T.
$$

These fronts form a **periodic ribbon of homotopic loops**.

A key property of stationary reentry is that each loop in the ribbon satisfies

$$
L_g(\Gamma_t) = T.
$$

Therefore individual wavefronts do not correspond to isolated critical points of the length functional.

Instead, they should be viewed as a **foliation of loops organized around a central geometric cycle**.

## 5. Organizing Cycles and Wave Ribbon Structure

We introduce the concept of an **organizing cycle**.

An organizing cycle is a closed loop

$$
\gamma : S^1 \to M
$$

around which a periodic pattern of activation wavefronts is organized.

In a stable reentrant regime, the observed wavefronts do not coincide with a single fixed loop. Instead, they form a **family of curves evolving in time**, which can be viewed as a ribbon surrounding the organizing cycle.

### 5.1 Wave Ribbon Structure

Let

$$
\Gamma_t \subset M
$$

denote the activation front at time $t$. In a periodic reentrant regime with period $T$, these fronts satisfy

$$
\Gamma_{t+T} = \Gamma_t.
$$

The set of curves $\{\Gamma_t\}_{t\in[0,T]}$ forms a **wave ribbon** surrounding the organizing cycle. Each $\Gamma_t$ is a closed curve homotopic to the organizing cycle $\gamma$.

Geometrically, this structure can be interpreted as a **foliation of a tubular neighborhood of&#x20;****$\gamma$****&#x20;by closed loops**. The organizing cycle therefore represents the central trajectory around which the wavefronts circulate.

This ribbon structure explains why individual wavefronts should not be identified with critical points of the length functional: the periodic solution corresponds to a **family of loops with identical travel time**, rather than to a single distinguished loop.

### 5.2 Examples

Typical examples of organizing cycles include:

* anatomical macro-reentry around scars or valves,
* the core circle of a solid torus,
* the trajectory of the tip of a spiral wave.

In each case, the wavefronts propagate as a ribbon surrounding a central geometric trajectory.

The organizing cycle therefore acts as the **geometric backbone of the periodic activation pattern**, capturing the topology of the reentrant circuit independently of the detailed shape of individual wavefronts.

## 6. Variational Motivation from the Eikonal Regime

The introduction of a variational description of organizing cycles is motivated by the **eikonal approximation** of wave propagation.

In this regime, activation times satisfy an anisotropic eikonal equation of the form

$$
g^{-1}(d\tau,d\tau)=1,
$$

where $\tau(x)$ denotes the activation time and $g$ is the electrophysiological metric introduced in Section 2.

This formulation identifies the metric $g$ as the geometric structure governing propagation. In particular, locally preferred propagation trajectories correspond to **geodesics of the metric&#x20;****$g$**.

For a stable periodic reentry, the organizing cycle must be globally self-consistent from one period to the next. This suggests a closed-loop analogue of **Fermat’s principle**: the organizing cycle should be stationary for the total electrophysiological travel time within its free homotopy class.

This naturally leads to the travel-time functional

$$
L_g(\gamma)
=
\int_0^1 \|\dot{\gamma}(t)\|_g \, dt .
$$

In practice it is convenient to work with the associated energy functional

$$
E(\gamma)
=
\frac12
\int_0^1
\|\dot{\gamma}(t)\|_g^2 dt .
$$

Although the full electrophysiological dynamics is not a gradient flow of $E$, the eikonal regime provides a natural variational motivation for modeling organizing cycles as **critical loops of the travel-time functional&#x20;****$L_g$**, or equivalently of $E$ after fixing parametrization.

## 7. Variational Formulation

We therefore consider the energy functional

$$
E : \Lambda M \to \mathbb{R}
$$

defined by

$$
E(\gamma)
=
\frac12
\int_0^1
\|\dot{\gamma}(t)\|_g^2 dt .
$$

Here

$$
\Lambda M = H^1(S^1,M)
$$

denotes the free loop space of $M$.

The functional $E$ is $C^2$ on the Hilbert manifold $\Lambda M$.

Critical points of $E$ satisfy the geodesic equation

$$
\nabla_t \dot{\gamma} = 0.
$$

Thus critical loops of $E$ correspond to **closed geodesics of the electrophysiological metric&#x20;****$g$**.

Within the present framework, organizing cycles of stable reentrant patterns are modeled as such closed geodesics.

## 8. Morse–Bott Structure of the Energy Functional

In classical Morse theory, critical points are typically assumed to be non-degenerate and isolated. However, the energy functional $E$ on the free loop space $\Lambda M$ possesses a natural **$S^1$****-symmetry** arising from time translation along loops: if $\gamma(t)$ is a closed geodesic, then $\gamma(t+\theta)$ is also a closed geodesic for any $\theta \in S^1$.

Consequently, critical points of $E$ are not isolated but appear as **critical manifolds**, typically $S^1$-orbits of closed geodesics. This naturally leads to a **Morse–Bott framework** rather than classical Morse theory.

### 8.1 Definition of Critical Manifolds

A connected submanifold $C \subset \Lambda M$ is called a **Morse–Bott critical manifold** of $E$ if:

1. Every point $\gamma \in C$ is a critical point of $E$ (i.e., a closed geodesic).
2. The Hessian $d^2E_\gamma$ is non-degenerate in directions normal to $C$.

In the context of cardiac reentry, each such manifold $C$ corresponds to a **family of organizing cycles** related by temporal phase shifts. The $S^1$ orbit within $C$ represents the temporal translation symmetry of the circulating wave along its organizing cycle.

### 8.2 Index and Physical Interpretation

For a Morse–Bott critical manifold $C$, the **Morse index** $\lambda(C)$ is defined as the dimension of the maximal subspace of the normal bundle on which the Hessian $d^2E$ is negative-definite.

Although the electrophysiological dynamics is not strictly a gradient flow of $E$, the Morse index provides a useful heuristic interpretation:

* **$\lambda(C)=0$****&#x20;(local minima):** organizing cycles that are expected to be dynamically stable under small perturbations.
* **$\lambda(C)>0$****&#x20;(saddle-type manifolds):** unstable periodic configurations that may act as transition states between different reentrant regimes.

### 8.3 Topological Constraints

Morse–Bott inequalities relate the topology of the loop space $\Lambda M$ to the indices of the critical manifolds:

$$
\sum_C t^{\lambda(C)} P_t(C) \ge P_t(\Lambda M)
$$

where the sum runs over all Morse–Bott critical manifolds of $E$.

These relations suggest that the **topological complexity of the myocardial domain** influences the number of organizing cycles. In particular, the non-trivial topology of $\Lambda M$, induced by the fundamental group $\pi_1(M)$, leads variationally to closed geodesic representatives in many homotopy classes.

## 9. Minimal Organizing Cycles

For each free homotopy class $\alpha \in [S^1,M]$, we consider the minimal electrophysiological length

$$
\ell(\alpha)
=
\inf_{\gamma \in \Lambda_\alpha M} L_g(\gamma).
$$

This quantity represents the shortest possible conduction time for an organizing cycle within the topological class $\alpha$.

### 9.1 Existence of Minimal Representatives

Under mild regularity assumptions on the domain $(M,g)$, variational arguments imply the existence of a loop

$$
\gamma_\alpha \in \Lambda_\alpha M
$$

such that

$$
L_g(\gamma_\alpha) = \ell(\alpha).
$$

In particular, if $M$ is compact and the boundary $\partial M$ is non-penetrable for wave propagation (for instance representing anatomical barriers or ablation lines), minimizing sequences remain confined in the interior of $M$.

In this case, standard compactness arguments in the Sobolev loop space $H^1(S^1,M)$ yield the existence of a minimizing loop.

### 9.2 Geodesic Character of Minimal Cycles

Any minimizing loop $\gamma_\alpha$ is stationary for the length functional, and therefore satisfies the geodesic equation

$$
\nabla_t \dot{\gamma}_\alpha = 0.
$$

Thus minimal organizing cycles correspond to **closed geodesics of the electrophysiological metric**.

### 9.3 Interpretation

The quantity

$$
\ell(\alpha)
$$

can therefore be interpreted as the **minimal organizing period** associated with the topological class $\alpha$.

If several classes satisfy the refractory compatibility condition

$$
\ell(\alpha) > R_{\mathrm{eff}},
$$

then multiple geometrically distinct organizing cycles may coexist in the domain.

## 10. Refractory Constraint and Admissible Cycles

In addition to geometric constraints, the existence of a sustained reentrant circuit requires compatibility with the refractory properties of the tissue.

Let

$$
R_{\mathrm{eff}}(x,T)
$$

denote the **effective refractory time** at position $x \in M$, which may depend on the local activation period $T$ due to restitution effects.

For a candidate organizing cycle $\gamma$, the total conduction time is

$$
T(\gamma) = L_g(\gamma).
$$

For the reentrant wave to propagate continuously, each point of the cycle must have recovered excitability when the wave returns. This leads to the **refractory compatibility condition**

$$
T(\gamma) > R_{\mathrm{eff}}(x,T(\gamma))
\quad \text{for all } x \in \gamma.
$$

A sufficient condition for admissibility is therefore

$$
L_g(\gamma) > \sup_{x \in \gamma} R_{\mathrm{eff}}(x,L_g(\gamma)).
$$

### 10.1 Minimal Organizing Length

For each free homotopy class $\alpha$, define

$$
\ell(\alpha)
=
\inf_{\gamma \in \Lambda_\alpha M} L_g(\gamma).
$$

A necessary condition for the existence of a reentrant circuit in class $\alpha$ is

$$
\ell(\alpha) > R_{\mathrm{crit}},
$$

where

$$
R_{\mathrm{crit}} = \sup_{x \in M} R_{\mathrm{eff}}(x,\ell(\alpha)).
$$

Classes satisfying this inequality constitute the set of **refractory-admissible organizing cycles**.

### 10.2 Physiological Interpretation

This constraint corresponds to the classical requirement that the **cycle length of the reentrant circuit must exceed the effective refractory period** along the entire pathway.

In electrophysiological terms, the condition expresses the balance between:

* the **geometric conduction time** determined by the electrophysiological metric $g$, and
* the **local recovery dynamics** of the tissue.

Stable reentry emerges only when the geometric cycle length is sufficiently long to allow full recovery before the next activation.

## 11. Examples: 

### 11.1 Solid Torus

Let

$$
M = D^2 \times S^1.
$$

Then

$$
\pi_1(M) \cong \mathbb{Z}.
$$

Thus there is a single primitive class of loops.

In an axisymmetric geometry the **core circle**

$$
\gamma_0
$$

is the unique primitive closed geodesic in the interior.

Its length

$$
L_0 = L_g(\gamma_0)
$$

determines the organizing cycle of the reentrant ribbon.

Wavefronts form a periodic ribbon surrounding this core circle.



### 11.2 Genus-2 Handlebody

Consider a genus-2 handlebody

$$
H_2.
$$

Its fundamental group is

$$
\pi_1(H_2) \cong F_2,
$$

the free group on two generators $a$ and $b$.

Possible organizing cycles correspond to classes

$$
a,\quad b,\quad ab,\quad ab^{-1},\dots
$$

Each class has an associated minimal organizing length

$$
\ell(a),\quad \ell(b),\quad \ell(ab),\dots
$$

In symmetric geometries one may have

$$
\ell(a) = \ell(b),
$$

leading to **multiple competing organizing cycles**.

## 12. Stability and Transition Energy

The stability of a cycle $\gamma$ is determined by the second variation of the energy:

$$
\delta^2 E_\gamma .
$$

If the second variation is positive in transverse directions, the cycle is locally stable.

Transitions between circuits correspond to trajectories crossing saddle-type critical manifolds.

The associated energy barrier

$$
\Delta E
$$

can be interpreted as a **transition energy** between organizing cycles.

## 13. Effect of Ablation

Within the present framework, ablation can be interpreted geometrically as a modification of the electrophysiological domain $(M,g)$.

Two distinct mechanisms arise.

### 13.1 Topological Modification

Ablation lesions may alter the topology of the excitable domain,

$$
M \rightarrow M'.
$$

This modification changes the set of free homotopy classes $[S^1,M]$.

In particular, some classes of organizing cycles may disappear entirely. For example, a linear lesion connecting two anatomical boundaries may eliminate a previously existing reentrant circuit by removing its homotopy class.

### 13.2 Modification of the Electrophysiological Metric

Ablation also alters local conduction properties, which can be modeled as a perturbation of the electrophysiological metric

$$
g \rightarrow g'.
$$

This modification changes the electrophysiological length

$$
L_g(\gamma)
$$

and therefore the minimal organizing lengths

$$
\ell(\alpha).
$$

As a consequence, previously admissible organizing cycles may become incompatible with the refractory constraint

$$
\ell(\alpha) > R_{\mathrm{eff}}.
$$

### 13.3 Variational Landscape

Since the energy functional

$$
E(\gamma)
=
\frac12 \int_0^1 \|\dot{\gamma}(t)\|_g^2 dt
$$

depends on both the geometry of $M$ and the metric $g$, ablation effectively modifies the **variational landscape** of organizing cycles.

Possible effects include:

* disappearance of a stable organizing cycle,
* emergence of new competing cycles,
* modification of transition barriers between circuits.

  In this sense, ablation may be interpreted as a **geometric surgery on the space of organizing cycles**.

## 14. Topological Constraints

The topology of the loop space $\Lambda M$ constrains the structure of the energy functional.

In particular, Morse–Bott theory implies that the number and type of organizing cycles are influenced by the topology of the underlying domain.

This provides a geometric explanation for the multiplicity of possible reentrant circuits in complex anatomical substrates.

## 15. Discussion

The central conceptual point of this framework is the distinction between

* **wavefronts**, which form periodic ribbons of loops, and
* **organizing cycles**, which act as geometric backbones of the dynamics.

The latter can be naturally modeled as closed geodesics of the electrophysiological metric.

This viewpoint provides a unified geometric interpretation of

* circuit selection
* multistability of reentry
* transitions between circuits
* effects of ablation.











# Potential directions (to be discussed)

## Refractory Constraint as an Effective Variational Penalty

*This may be interesting compared to a standard "static" model, namely to account for scar?*

The simple condition

$$
L_g(\gamma) > R_{\mathrm{eff}}
$$

provides a useful first approximation for reentry admissibility, but it treats refractoriness as a hard threshold. A more physiological description is obtained by incorporating recovery directly into the variational framework.

For a candidate organizing cycle $\gamma$, let

$$
T(\gamma) = L_g(\gamma)
$$

denote its electrophysiological period. Define the local recovery margin

$$
m_\gamma(x)
=
T(\gamma) - R_{\mathrm{eff}}(x,T(\gamma)).
$$

Reentry is locally admissible where $m_\gamma(x)>0$.

We then define a refractory penalty functional

$$
\mathcal{R}(\gamma)
=
\int_\gamma
\bigl[
R_{\mathrm{eff}}(x,T(\gamma)) - T(\gamma)
\bigr]_+
\, ds_g,
$$

where $[u]_+ = \max(u,0)$ denotes the positive part.

This term vanishes when the cycle is fully compatible with local recovery, and becomes positive when part of the cycle remains refractory at the time of return.

An effective functional for organizing cycles is therefore

$$
\mathcal{F}(\gamma)
=
L_g(\gamma)
+
\mu\,\mathcal{R}(\gamma),
$$

with $\mu>0$ a penalty parameter.

In this formulation, admissible organizing cycles are modeled as local minimizers (or more generally critical manifolds) of $\mathcal{F}$ rather than of $L_g$ alone. This yields a more natural geometric-physiological compromise: stable reentry selects cycles that are both short in electrophysiological length and compatible with tissue recovery.

