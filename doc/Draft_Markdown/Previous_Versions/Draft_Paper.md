---
title: "An analytical approach to transient dynamics in predator-prey systems"
subtitle: ""
author: Sergio Cobo-Lopez, Matthew Witt, Forest L. Rohwer, and Toni Luque
csl: "podzemna-voda.csl"
abstract: "This paper presents a mathematical tool to predict tipping points in transient dynamical systems. Transient dynamics are characterized by shifts between two consecutive quasi-stable dynamics. These shifts often occur very fast compared to the system's characteristic timescale and cannot be predicted by simply looking at the dynamics. Our method stems from the assumption that given a system of differential equations describing a dynamical system, not all the terms in the equations are active all the time. Terms are activated and inactivated by tipping points: critical values of the state variables in the system. As a proof of concept, we apply our method to a predator-prey model of bacteria and bacteriophage, the most abundant biological agents on earth. In this case, tipping points activating or inactivating ecological terms are given by critical concentrations of phage and bacteria. Our method identifies the conditions for transient dynamics and successfully predicts the tipping points leading to bacterial or phage extinction. Our results show that different characteristic timescales of bacteria and phage result in transient dynamics."
graphics: yes
header-includes:
  - \usepackage{color}
  - \usepackage{caption}
  - \usepackage{anysize}
  - \usepackage{amsmath}
  - \usepackage{siunitx}

---

[comment]: <> (To compile this document with a bibliography: pandoc Draft_Paper.md -o Draft_Paper.pdf --bibliography Transient_Dynamics.bib)

## Introduction

### Hypothesis and background: All systems are transient

"In all of my universe I have seen no law of nature, unchanging and inexorable. This universe presents only changing relationships which are sometimes seen as laws by short-lived awareness. [...]. If you must label the absolute, use its proper name: Temporary" (page 565, The stolen journals section, God Emperor of Dune by Frank Herbert). All dynamical systems in nature display transient dynamics. Systems that appear to behave asymptotically are simply transient over long timescales. The Solar System is an excellent example of this behavior: although it is predicted to be stable on a scale of a couple billion years [@Hayes2007], it has been proved to have a chaotic behavior on longer timescales. For instance, [@Laskar2009] found a small probability of a collision between the Earth and Mars, Venus, or Mercury over the next 3.4 billion years. In general, during the occurrence of transient dynamics, the variables of the system can show a gradual change or a minor drift followed by a sudden and extensive change [@rocha_peterson_bodin_levin_2018]. Qualitatively, these transient regimes can be due to the competition between underlying equilibria [@ludwig_jones_holling_1978], [@doi:10.1126/science.aat6412] to changes in the physical parameters in the system [@ludwig_jones_holling_1978], [@doi:10.1126/science.aat6412], or to competing time scales associated to different agents in the system [@ludwig_jones_holling_1978], [@cairns:2009],[@ROACH201738]. Unlike the study of asymptotic dynamics or equilibria, the quantitative framework to classify and investigate transient dynamics is still under development [@9910223127702121] [@doi:10.1126/science.aat6412].

### What has been done so far

A common approach to study transient dynamics is to assume that the system is stable in the vicinity of the equilibrium. However, this is not always true. Another practical approach is to approximate a physical system using sets of ordinary differential equations and integrate them numerically to obtain the trajectories in the regions of interest [@bashkirtseva_ryashko_2018], [@cairns:2009], [@cushing_dennis_desharnais_costantino_1998], [@ROACH201738], [@rabinovich_varona_selverston_abarbanel_2006], [@rinaldi_scheffer_2000]. Fitting these models to empirical data provides parameter values that can reproduce the transient dynamics and can be used to interpret experiments and make predictions. This approach, however, requires an exhaustive exploration of the parameter space to identify alternative situations associated to different regime shifts, which becomes particularly challenging as additional mechanisms and dynamical variables are included to reflect the complexity of the system of study [@cairns:2009].

### How to solve the problem

Here, we present a framework to study transient dynamics that is agnostic about equilibrium conditions, and eases the problem of model complexity. The most important assumption in our approach is that, in general, not all the mechanisms in a dynamical system are active all the time, which reduces the model complexity. Our approach uses the activation/inactivation of mechanisms to predict transient dynamics. We use per capita rates [@49249261], to estimate when a mechanism is inactive or not. In this paper, we apply this approach to a predator-prey system that models the ecology of a population of *E. coli* bacteria and virulent phage T4. We chose this system because transient dynamics are pervasive in ecology (extinctions or migrations are typical examples thereof), and bacteria and bacteriophages are the most abundant biological entities on Earth. Specifically, the combination of T4 phage and *E. Coli* is very common in experiments that can benefit from accurate mathematical models

"Everything flows. In this astonoshing and unfamiliar landscape, everything -the snap of a grasshopper, a rustling cottonwood, an owl hooting at night- brings us into the moment and pulls us into the flow" Heraclitus

"In all of my universe I have see no law of nature, unchanging and inexorable. This universe presents only changing relationships which are sometimes seen as laws by short-lived awareness. [...]. If you must label the absolute, use its proper name: Temporary" page 565, The stolen journals section, God Emperor of Dune by Frank Herbert).

All dynamical systems in nature display transient dynamics. Systems that appear to be stable or in the equilibrium are simply transient over long timescales. The Solar System is an excellent example of this behavior: although it is stable on a scale of a couple billion years [@Hayes2007], it has been proved to have a  chaotic behavior on longer timescales. For instance, [@Laskar2009] found a small probability of a collision between the Earth and Mars, Venus, or Mercury over the next 3.4 billion years. In general, during the occurrence of transient dynamics, the variables of the system can show a gradual change or a minor drift followed by a sudden and extensive change [@rocha_peterson_bodin_levin_2018]. Qualitatively, these transient regimes can be due to the competition between underlying equilibria [@ludwig_jones_holling_1978], [@doi:10.1126/science.aat6412] to changes in the physical parameters in the system [@ludwig_jones_holling_1978], [@doi:10.1126/science.aat6412], or to competing time scales associated to different agents in the system [@ludwig_jones_holling_1978], [@cairns:2009],[@ROACH201738]. Unlike the study of asymptotic dynamics or equilibria, the quantitative framework to classify and investigate transient dynamics is still under development [@9910223127702121] [@doi:10.1126/science.aat6412].

### What has been done so far

A common approach to study transient dynamics is to assume that the system is stable in the vicinity of the equilibrium. As explained above, this not need to be the case. Another practical approach is to approximate a physical system using sets of ordinary differential equations and integrate them numerically to obtain the trajectories in the regions of interest [@bashkirtseva_ryashko_2018], [@cairns:2009], [@cushing_dennis_desharnais_costantino_1998], [@ROACH201738], [@rabinovich_varona_selverston_abarbanel_2006], [@rinaldi_scheffer_2000]. Fitting these models to empirical data provides parameter values that can reproduce the transient dynamics and can be used to interpret experiments and make predictions. This approach, however, requires an exhaustive exploration of the parameter space to identify alternative situations associated to different regime shifts, which becomes particularly challenging as additional mechanisms and dynamical variables are included to reflect the complexity of the system of study [@cairns:2009].

### How to solve the problem

Here, we present a framework for the study of transient dynamics that is agnostic about equilibrium conditions, and eases the problem of model complexity. The most important assumption in our approach is that, in general, not all the mechanisms in a dynamical system are active all the time, which reduces the model complexity. Our approach uses the activation/inactivation of mechanisms to predict transient dynamics. We use per capita rates [@49249261], to estimate when a mechanism is inactive or not. In this paper, we apply this approach to a predator-prey system that models the ecology of a population of *E. coli* bacteria and virulent phage T4. We chose this system because transient dynamics are pervasive in ecology (extinctions or migrations are typical examples thereof), and bacteria and bacteriophages are the most abundant biological entities on Earth. Specifically, the combination of T4 phage and *E. Coli* is very common in experiments that can benefit from accurate mathematical models


## Results

* We model the ecology of bacteria and virulent phages with Lotka-Volterra equations. The equations encode four mechanistic terms: bacterial growth, phage predation, viral burst, and viral decay (Figure 1 a, Equation 1). We chose parameters that correspond to an *E. coli* bacteria and a T4 virulent phage (Table 1), a very common experimental setup. The dynamics of this system are characterized by a very fast bacterial growth and a slow phage decay.

* We ran the model for 14 hours, a typical time in experiments. The result of the simulation (Figure 1 b) shows the dynamics of *E. coli* (blue) and T4 (red) in a logarithmic scale. The concentration of *E. coli* grows exponentially over the first eight hours, and drops to extinction afterward. T4 remain constant during the first four hours, then grow as a double exponential for another four hours, and it plateaus afterward.

* To get a better insight of these dynamics, we look at the contribution of each term separately (Figure 1 c). The contributions are normalized by the total concentration of bacteria and phage and referenced to the fastest timescale of the system. In this case, the characteristic timescale is given by the bacterial growth rate $r$:
=======
[comment]: <> (*Some ecologists have proposed theoretical frameworks that can be used to study transient phenomena in population dynamics. Specifically, [@49249261] discusses population dynamics from first principles and argues that ecological processes are driven by individual processes such as births, predations, or deaths. [@49249261] looks at population dynamics from a per capita perspective: how many individuals are added or removed from the system per individual. Although per capita rates reference individual changes to collective changes, they do not account for the different timescales that operate in an ecosystem.)

* We used Lotka-Volterra equations to model the ecology of bacteria and virulent phages. The equations encode four mechanistic terms: bacterial growth, phage predation, viral burst, and viral decay (Figure 1 a, Equation 1). We chose parameters that correspond to an *E. coli* bacteria and a T4 virulent phage (Table 1), a very common experimental setup. The dynamics of this system are characterized by a very fast bacterial growth and a slow phage decay.

* We ran the model for 14 hours, a typical time in experiments (Figure 1 b). The dynamics of this system (Figure 1b) show a steady growth of *E. coli* followed by a rapid extinction after eight or nine hours. The concentration of T4 remains almost constant for the first four or five hours of the experiment, then grows very fast for three hours, and saturates afterwards.

* To get a better insight of these dynamics, we look at the contribution of mechanistic terms separately (Figure 1 c). We normalized these terms by the total concentration of bacteria and phage (see Equation 1) and obtain per capita contributions. Then we referenced per capita contributions to a timescale. In this case, we set the bacterial growth rate $r$ as the reference timescale, because it is the fastes process in the system: 

\begin{align}
\frac{dB}{dt}&= rB - d P B & \frac{1}{B}\frac{dB}{dt}&= r - d P &  \frac{1}{rB}\frac{dB}{dt}&= 1 - \frac{d}{r} P  \\
\frac{dP}{dt}&= cdBP - m P & \frac{1}{P}\frac{dP}{dt}&= cdB - m       &  \frac{1}{rP}\frac{dP}{dt}&= \frac{cd}{r}B - \frac{m}{r} \nonumber
\end{align}

* Introduce and justify tipping points or critical concentrations.

These time aggregated per capita rates give the number of bacteria or phage that enter or leave the system per unit of bacteria or phage over a time frame equivalent to the duplication time of the average bacteria. Figure 1 C shows these rates over time. While the contributions of bacterial growth and decay are constant, bacterial predation and viral burst change over time and, particularly, change with phage and bacterial density, respectively. More specifically, we observe major changes in these rates around the point where bacteria start decreasing. From this observations, we conclude that mechanistic terms have different weights at different times. 

\begin{figure}[htbp]\centering\includegraphics[width=0.5\textwidth, height=!]{/home/sergio/work/Github/needle-finder/docs/Finished_Figures/Paper_Figures/Figure_1.pdf}\caption{ \textbf{Identification of active terms in Lotka-Volterra dynamics with aggregated per capita rates}. \textbf{a)} Graphical representation of the predator-prey model in a phage bacterial system. \textbf{b)} Dynamics of the predator and prey concentrations when r is the dominant timescale. \textbf{c)} Contribution of terms to rates. Aggregated per capita rates are obtained by normalizing rates by the total population of prey or predator and the dominant timescale (r). Deeper red and blue represent increasing contributions of predator and prey-like terms to the dynamics, respectively.}\end{figure}

* We hypothesize that we can inactivate certain terms in our model if their contribution is below a certain critical level and still explain the dynamics of the system accurately. Therefore we can build a minimal model: a concatenation of simplified dynamics or versions of the Lotka-Volterra equations with inactive or negiglible processes excluded. If the aggregated per capita contribution of a term is smaller than a critical value $\epsilon$, then that term is inactive.

* Without loss of generality, let us assume that $\epsilon=1$. That means that any mechanistic process that introduces or removes less than one bacteria or phage per capita over the relevant timescale will be inactive.

* We can use the predictive power of time aggregated contributions per capita to predict transient dynamics. In particular, we can build a minimal model: a concatenation of simplified dynamics or versions of the Lotka-Volterra equations with inactive or negiglible processes excluded. We use a parameter $\epsilon$ to control the critical contribution that sets a term inactive. Figure 2 a shows the simplified dynamics of the system described above when $\epsilon=1$: first, the only term active is growth (of bacteria). When the bacteria reach the critical concentration, viral burst becomes active. Viral burst increases phage concentration, and when phage reach the critical concentration, they activate the predation. As long as the burst is active, predation will keep increasing and it eventually outgrowths bacterial growth. That induces a decrease of bacteria that will inactivate the burst.

* Figure 2 b shows the results of the full model and the minimal model for $\epsilon=1$. The star, circle, and triangle, indicate the times at which the critical concentrations $B_c$ and $P_c$ were reached. Roman numbers indicate the simplified dynamics as showed in Figure 2 a. The minimal model reproduces the full dynamics very well. It understimates viral growth at the end of Simplified dynamics I because it is not considering the viral burst.

* Firgure 2 c shows the number of mechanistic processes active over time. Importantly, there is no point at which all four processes are active. The viral decay never activates, because the it is very slow compared to the bacterial growth.

\begin{figure}[htbp]\centering\includegraphics[width=0.5\textwidth, height=!]{/home/sergio/work/Github/needle-finder/docs/Finished_Figures/Paper_Figures/Figure_2.pdf}\caption{ \textbf{Minimal model for r dominant timescale}. \textbf{a)} Simplified dynamics activate and inactivate terms as threshold concentrations of predator and prey are crossed.  \textbf{b)} Full dynamics (red and blue solid lines) versus minimal model (dashed lines). Dotted vertical lines indicate critical concentrations that drive the change in simplified dynamics depicted in \textbf{a)}. \textbf{c)}, Number of active terms over time.}\end{figure}


* We next analyze the case where the viral decay rate is the dominant timescale. Although not common, it is theoretically possible and worth analyzing from our methodology. Figure 3 a shows the simplified dynamics in this case. 

\begin{figure}[htbp]\centering\includegraphics[width=0.5\textwidth, height=!]{/home/sergio/work/Github/needle-finder/docs/Finished_Figures/Paper_Figures/Figure_3.pdf}\caption{ \textbf{Minimal model for m dominant timescale}. \textbf{a)} Sequence of simplified dynamics. \textbf{b)} Full dynamics (red and blue solid lines) versus minimal model (dashed lines). Dotted vertical lines indicate critical concentrations that drive the change in simplified dynamics depicted in \textbf{a)}. \textbf{c)}, Number of active terms over time. }\end{figure}

Finally, we analyze the case where timescales are balanced and $$r=m$$

\begin{figure}[htbp]\centering\includegraphics[width=0.75\textwidth, height=!]{/home/sergio/work/Github/needle-finder/docs/Finished_Figures/Paper_Figures/Figure_4.pdf}\caption{ \textbf{Minimal model for m dominant timescale}. \textbf{a)} Sequence of simplified dynamics for $\epsilon=1$ ($\epsilon=0.1$ is confined to dynamic I). \textbf{b)} and \textbf{d)}  Full dynamics (red and blue solid lines) versus minimal model (dashed lines) for $\epsilon=1$ and $\epsilon=0.1$. Dotted vertical lines indicate tipping points between different dynamics. \textbf{c)} and \textbf{e)}, Number of active terms over time.}\end{figure}

\begin{figure}[htbp]\centering\includegraphics[width=0.75\textwidth, height=!]{/home/sergio/work/Github/needle-finder/docs/Finished_Figures/Paper_Figures/Figure_4_Draft_1.pdf}\caption{ \textbf{Minimal model for m dominant timescale}. \textbf{a)} Sequence of simplified dynamics for $\epsilon=1$ ($\epsilon=0.1$ is confined to dynamic I). \textbf{b)} and \textbf{d)}  Full dynamics (red and blue solid lines) versus minimal model (dashed lines) for $\epsilon=1$ and $\epsilon=0.1$. Dotted vertical lines indicate tipping points between different dynamics. \textbf{c)} and \textbf{e)}, Number of active terms over time.}\end{figure}


## Discussion


## Methods

Transient dynamics are pervasive in ecology: extinctions or migrations are typical examples of ecological shifts that dramatically change the dynamics of ecosystems. Some ecologists have proposed theoretical frameworks that can be used to study transient phenomena in population dynamics. Specifically, [@49249261] discusses population dynamics from first principles and argues that ecological processes are driven by individual processes such as births, predations, or deaths. [@49249261] looks at population dynamics from a per capita perspective: how many individuals are added or removed from the system per individual. Although per capita rates reference individual changes to collective changes, they do not account for the different timescales that operate in an ecosystem. Here, we develop a method that combines per capita rates and relevant timescales to predict transient dynamics in a predator-prey model (Lotka-Volterra equations). Specifically, we apply our method to the interaction of bacteria and bacteriophages, Earth's two most abundant biological agents.             

In general, per capita rates are obtained by normalizing Lotka-Volterra equations by the total concentration of prey (bacteria) and predator (phage) (see Eq. 1). Our method goes a step further and normalizes per capita concentrations by the relevant timescale, the typical time of the fastest ecological process. In our case, the relevant timescales are determined either by the bacterial growth rate $r$ or the phage decay rate $m$. This normalization gives a time aggregated per capita rate that represents the number of bacteria/phage that are added or removed from the system per unit of bacteria/phage over the relevant timescale. Looking at the ecosystem from this perspective, we can track the impact of individual rates on collective dynamics and predict transient dynamics. More specifically, we define critical concentrations of bacteria or phage that activate or inactivate mechanistic terms in the Lotka-Volterra equations.

Our results suggest that a clearly dominant timescale leads to transient dynamics characterized by phage or bacterial extinction and that quasi-stable dynamics are only going to be observed when timescales are balanced. It is important to recall that, by our own assumption, quasi-stability can only mean that the model is incomplete: because all dynamical systems are transient, there might be a missing term that accounts for transients, even if that term only operates over long timescales.

### How did I do it

Without loss of generality, we chose parameter values for the Lotka-Volterra equations that simulate a system of *E. coli* bacteria and virulent T4 phage (see Table 1). These parameters account for the four mechanistic terms in the model: bacterial growth rate, phage predation, viral burst, and phage decay (see Figure 1,(a).) We took two different values for $r$ and $m$ because we wanted to simulate three different scenarios: two in which there is a dominant timescale ($r>m$ or $m>r$) and one in which timescales are equivalent ($r=m$).

| Parameter | Description | Values | Source|
| ------- | ----------------- | ------------------ | ---------------- | 
| r | Maximum growth rate |0.9 \si{h^{-1}}, $\num{3.7e-5} \si{h^{-1}}$ | [@https://doi.org/10.1111/1462-2920.15640] |
| d | Infection rate |\num{3e-8} \si[per-mode=symbol]{\ml\per\hour}|[@10.1371/journal.pbio.0040193]|
| c | Burst size | 150 |[@10.1371/journal.pbio.0040193]  |
| m | Decay rate |\num{2.8e-3} \si{h^{-1}}, 0.1 \si{h^{-1}} | [@10.1371/journal.pbio.0040193], [@doi:10.1128/aem.58.11.3721-3729.1992]  |
Table: 1

Let us consider first the situation in which $r$ is the dominant timescale. Then, we take per capita rates and normalize them by $r$ to obtain time aggregated per capita rates:

\begin{align}
\frac{dB}{dt}&= rB - d P B & \frac{1}{B}\frac{dB}{dt}&= r - d P &  \frac{1}{rB}\frac{dB}{dt}&= 1 - \frac{d}{r} P  \\
\frac{dP}{dt}&= cdBP - m P & \frac{1}{P}\frac{dP}{dt}&= cdB - m       &  \frac{1}{rP}\frac{dP}{dt}&= \frac{cd}{r}B - \frac{m}{r} \nonumber
\end{align}

Time aggregated per capita rates are just a different way to look at the Lotka-Volterra equations. These rates allow us to measure the impact of individual processes (over a specific time) on the bacterial and phage communities. Specifically, we can set a critical value $\epsilon$ for these rates, below which the corresponding mechanistic term is inactivated. Suppose, for instance, that $\epsilon=1$; in our case, that means that growth is always active because it is constant and equal to $1$. Decay is always inactive because $\frac{m}{r}=\num{3.11e-2}$<1. Finally,  predation and burst  depend on the concentrations of phage and bacteria, respectively. Their critical concentrations are:

\begin{equation*}
P_c=\epsilon \frac{r}{d} \qquad B_c=\epsilon \frac{r}{cd} ,
\end{equation*} 

$\epsilon$ ,thefore, acts as a tipping point for the mechanistic terms in the Lotka-Volterra equations. In our method, tipping points for time aggregated per capita rates activate or inactivate terms in the Lotka-Volterra equations. We can therefore solve the dynamics of this system with a minimal model: a concatenation of simplified dynamics or versions of the Lotka-Volterra equations with inactive terms excluded. Lower values of $\epsilon$ will generally result in minimal models closer to the original ones (more precise). In fact, we do not need to solve the equations; just by knowing the active and inactive terms, we can determine when and how the transient dynamics will occur.

It is important to recall that our method allows us to predict the outcome of the dynamics without even having to solve the predator-prey equations: when only growth is active, it is straightforward that the bacterial concentration will eventually reach the critical concentration $B_c$ and activate the burst (region II). The activation of the burst causes the increase of predator concentration which will reach $P_c$, thus activating the predation. Once the predation is active (region III), the bacterial concentration starts decaying. Because predation is a density-dependent term of the phage (predators), and because the phage concentration is increasing via the burst, predation will eventually outweigh growth. The concentration of bacteria will fall below $B_c$, thus inactivating the burst. Only growth and predation are active at this point (region IV). Because the predator concentration is effectively constant (decay is inactive), the predation will outweigh growth, and bacteria will go extinct.

In addition to that, the method can produce a minimal model.

Something similar happens when $m$ is the dominant timescale ($m>r$, $r=\num{3.7e-5} \si{h^{-1}}$, and $m=0.1 \si{h^{-1}}$). In this case, we have:

\begin{align}
\frac{1}{mB}\frac{dB}{dt}&= \frac{r}{m} - \frac{d}{m} P \\
\frac{1}{mP} \frac{dP}{dt}&=\frac{cd}{m}B - 1 \quad, &  \nonumber
\end{align}

where decay and growth are now active and inactive, respectively.
