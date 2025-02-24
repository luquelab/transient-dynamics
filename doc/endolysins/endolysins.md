---
<<<<<<< HEAD
title: "A model for endolysins in transients dynamics"
subtitle: ""
author: "Sergio Cobo-Lopez"
csl: "podzemna-voda.csl"
abstract: "test"
=======
title: "A model for endolysins in transients dynamics
subtitle: ""
author: 
csl: "podzemna-voda.csl"
abstract: 
>>>>>>> 9944fffb8a9fa6512692733ae4d9ae95bfad9942
graphics: yes
header-includes:
  - \usepackage{color}
  - \usepackage{caption}
  - \usepackage{anysize}
  - \usepackage{amsmath}
  - \usepackage{siunitx}
<<<<<<< HEAD
---

[comment]: <> (To compile this document with a bibliography: pandoc endolysins.md -o endolysins.pdf --bibliography ../references/References.bib)
=======

---

[comment]: <> (To compile this document with a bibliography: pandoc Dendolysins.md -o endolysins.pdf --bibliography ../references/References.bib)
>>>>>>> 9944fffb8a9fa6512692733ae4d9ae95bfad9942

## The model
We propose a predator-prey model to describe the interaction of endolysins and bacteria. This model consists of six parameters and five mechanistic terms (see equation 1 below):

1. Bacterial logistic growth: bacteria duplicate until they reach the carrying capacity.
2. Bacterial predation by endolysins: n endolysins to find a bacterial cell and kill it. 
3. Endolysin flow: endolysins are introduced into the system at a certain rate.
4. Endolysin decay by predation: n endolysins leave the system whenever they bump into a bacterial cell.
4. Endolysin decay: endolysins will decay if they do not bump into their target bacteria.


\begin{align}
\frac{dB}{dt}&= \underbrace{r(1-B/K)B}_{\text{logistic growth}} - \underbrace{d B nE}_{\text{endolysin predation}} & \\
\frac{dE}{dt}&= \underbrace{fE}_{endolysin flow} - \underbrace{dBnE}_{\text{endolysin predation}} - \underbrace{m E}_{\text{endolysin decay}} & 
\end{align}



| Parameter | Units|Description | Values | Source|
| ------- | --------|------------ | ------------------ | ---------------- | 
| $r$ |$h^{-1}$ | Maximum growth rate |  | |
| $K$ |$cells/ml$ | Carrying capacity   | | |
| $d$ |$ml/h$| Infection rate      | | |
| $n$| |Number of endolysins to kill a bacteria | | |
| $f$ |$h^{-1}$ | Flow rate of endolysins | | |
| $m$ |$h^{-1}$| Decay rate of endolysins | | |
Table: 1
