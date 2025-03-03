---
layout: default
title: Approach
nav_order: 4
header-includes:
    - \usepackage{siunitx}
---

*"In all of my universe I have seen no law of nature, unchanging and inexorable. This universe 
presents only changing relationships which are sometimes seen as laws by short-lived awareness. [...]. If you must label the absolute, use its proper name: Temporary"*
(page 565, The stolen journals section, God Emperor of Dune by Frank Herbert).


## 1. Hyperion: First part

Our hypothesis is that, if we describe mathematically the mechanisms governing phage bacterial dynamics, we will observe that not all mechanisms are active all the time. The activation and inactivation of mechanisms generates changes in the dynamics via tipping points. Given this hypothesis, we can:

1. Predict tipping points
2. Prevent tipping points
3. Identify conditions of stability

As a proof of concept, we model the ecology of *E. Coli* bacteria and T4 virulent phages with Lotka-Volterra equations:

$$
\begin{align}
\frac{dB}{dt}&= rB - d P B & \\
\frac{dP}{dt}&= cdBP - m P &  \nonumber
\end{align}
$$

These equations model four mechanistic terms: bacterial growth, phage predation, viral burst, and viral decay. We chose parameters that correspond to an *E. coli* bacteria and a T4 virulent \
phage (Table 1), a very common experimental setup.

| Parameter | Description | Values | Source|
| ------- | ----------------- | ------------------ | ---------------|
| r | Maximum growth rate |$$0.9 h^{-1}$$, $$3.7e-5 {h^{-1}}$$ | [@https://doi.org/10.1111/1462-2920.15640] |
| d | Infection rate |$\num{3e-8} \si[per-mode=symbol]{\ml\per\hour}$|[@10.1371/journal.pbio.0040193]|
| c | Burst size | 150 | [@10.1371/journal.pbio.0040193]  |
| m | Decay rate |$\num{2.8e-3} \si{h^{-1}}$, $0.1 \si{h^{-1}}$ | [@10.1371/journal.pbio.0040193], [@doi:10.1128/aem.58.11.3721-3729.1992]|

Table: 1 

We know want to determine the tipping points that activate or inactivate terms. The first step is to calculate the per-capita rates of each term; that is, how many bacteria or phage enter or leave the system per unit of bacteria or phage:

$$
\begin{align}
\frac{1}{B}\frac{dB}{dt}&= r - d P &  \\
\frac{1}{P}\frac{dP}{dt}&= cdB - m&  \nonumber
\end{align}
$$

Next, we fix a relevant timescale or observation time $$\tau_o$$. The timescale determines the mechanisms that will be active:

$$
\begin{align}
\frac{\tau_o}{B}\frac{dB}{dt}&= \tau_o r - \tau_o d P &  \\
\frac{\tau_o}{P}\frac{dP}{dt}&= \tau_o cdB - \tau_o m&  \nonumber
\end{align}
$$

We consider a mechanistic term $$P_i$$ to be active or observable if $$\tau_o P_i \geq 1$$. There are two types of processes: constant (growth or phage decay) and time dependent (predation or burst). In the latter case, the time dependence arises from the dependence on bacterial or phage concentrations. If a process is time dependent, the characteristic timescale will also be time-dependent and the process will have critical concentrations.

Given a timescale $$\tau_o$$, the conditions for the four mechanistic terms to be active are:

$$
\begin{align}
\tau_o r \geq 1 &  \\
\tau_o d P \geq 1 &, P_c=\frac{1}{\tau_o d} \nonumber \\
\tau_o cdB \geq 1 &, B_c=\frac{1}{\tau_o cd} \nonumber \\
\tau_o m \geq 1 &  \nonumber
\end{align}
$$


## Application to spoT mutant. Initial Elements:

We build a dynamic mechanistic model that includes the main biological mechanisms in the system under study:

$$
\begin{align*} 
\frac{dL}{dt}&=\underbrace{r\bigg(1 - \frac{L}{K}\bigg)}_{\text{Logistic growth}} - \underbrace{\mu_i L}_{\text{Induction}} - \underbrace{\delta dT}_{\text{Infection}}& \\
\frac{dT}{dt}&=\underbrace{\mu_i L}_{\text{Induction}} - \underbrace{m T}_{\text{Decay}} - \underbrace{d TL}_{\text{Infection}}&
\end{align*} $$

This model assumes the existence of two bioagents (Lysogens and Temperate phages) that can undergo the following processes: Lysogen logistic growth, lysogen induction, lysogen reinfection, viral burst, and viral decay. As a first approximation, let us assume that lysogens are totally immune against infections ($$\delta=0$$). This yields the simplified model:

$$
\begin{align*} 
\frac{dL}{dt}&=\underbrace{r\bigg(1 - \frac{L}{K}\bigg)}_{\text{Logistic growth}} - \underbrace{\mu_i L}_{\text{Induction}}& \\
\frac{dT}{dt}&=\underbrace{\mu_i L}_{\text{Induction}} - \underbrace{m T}_{\text{Decay}}&
\end{align*}$$

The processes considered in the model are controlled by the parameters shown in this table:

| Parameter | Description | Value| Source|
| ----------- | ----------- | ----------- | ----------- |
| r | Maximum Growth Rate |0.471 (wild type), 0.382 (spot- mutant)  | Experiment data |
| K | Carrying capacity |8.2x10 7 (Wild type), 1x10 8 (mutant)  | Concentrations of bacteria at the end of the experiments without CX  |
| $$\mu_i$$ | Induction rate |spoT-: 8.481e-12, WT + CX: 5.793e-09, spoT- + CX: 1.814e-09, WT: 3.993e-11 | Experiment results  |
| $$\delta$$ | Probability of lysogen infection | 0 |  |
| d | Infection rate | 0 |  |
| c | Burst size | 125 |Da. Paepe et al, 2006  |
| m | Decay rate | 0.012h<sup>-1</sup>| Da Paepe et al, 2006  |

Our first step will be to try to fit the experimental data to this model.
