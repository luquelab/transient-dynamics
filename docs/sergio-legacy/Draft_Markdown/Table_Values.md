---
title: "Values for important parameters"
graphics: yes
header-includes:
  - \usepackage{color}                                                                              
  - \usepackage{caption}                                                                            
  - \usepackage{anysize}                                                                            
  - \usepackage{amsmath}                                                                            
  - \usepackage{siunitx}                                                                            
                                                                                                    
---

[comment]: <> (To compile this document: pandoc Table_Values.md -o Table_Values.pdf)



|Dominant timescale| r  | m | r and m  |
| ------------ |----------------- | ----------------- | ------------------ | 
| $B_c=1/(\tau c d)$   | $15873.016$ cells/ml    | $256410.26$ cells/ml  | $1111.11$ cells/ml |
| $P_c=1/(\tau d)$ | $2380952.38$ virions/ml | $38461538.46$ virions/ml | $166666.66$ virions/ml |
| $B_0$   | $1000$ cells/ml    | $66666666.67$ cells/ml  | $55555.55$ cells/ml |
| $P_0$   | $1e4$ virions/ml    | $1$ virions/ml  | $3333333.33$ cells/ml |
| $\tau$ | $14$ h | $260$ h | $200$ h|
| $w_{growth}=r \tau$ | $12.6$ | $2.6e-08$  | $20.0$ |
| $w_{decay}=m \tau$ | $0.0042$ | $26.0$ | $20.0$ |
| Transition times |$t^c_{1}=3.07$ h,$t^c_{2}=7.79$ h,$t^c_{3}=9.068$ h |$t^c_{1}=19.50$ h,$t^c_{2}=33.98$ h, $t^c_{3}=81.46$ h | |
| Error Bacteria* | 0.0952 | 0.0295  |0.0 |
| SD Bacteria | 0.1448  | 0.0166  | 0.0 |
| Error Phage | 0.0321  | 0.0173  | 0.0 |
| SD Phage | 0.0241  | 0.0076  | 0.0 |
| Error Total | 0.0636  |0.0234 | 0.0  |
| SD Total| 0.0683  | 0.0114 | 0.0 |


*How the error was calculated:


1. For each point (time point) $i$ take the concentration of bacteria and phage given by the full model ($B^{f}_{i}$ $P^{f}_{i}$) and by the simplified model ($B^{m}_{i}$ $P^{m}_{i}$)
2. Calculate the absolute value of the relative error of bacteria  $| \epsilon_B| (i)=|\frac{B^{f}_{i} -B^{m}_{i} }{B^{f}_{i}}|$, the absolute value of the relative error of phage $| \epsilon_P| (i)=|\frac{P^{f}_{i} -P^{m}_{i} }{P^{f}_{i}}|$, and their average: $\epsilon_{Tot}(i) = 0.5 (| \epsilon_B| (i) + | \epsilon_P| (i))$
3. Iterate over all points. This gives the vectors $|\epsilon_B|=(| \epsilon_B| (1), | \epsilon_B| (2), \dots ,| \epsilon_B| (n))$, $|\epsilon_P|$, and $|\epsilon_{Tot}|$ 
4. The error of bacteria, phage, and the total error are the means of the vectors mentioned above: $\overline{|\epsilon_B|}$, $\overline{|\epsilon_P|}$, $\overline{|\epsilon_{Tot}|}$

