---
layout: default
title: Background
nav_order: 3
---

# Background

Hyperion is a modeling tool for the optimization of experimental measurements. In particular, Hyperion uses mechanistic dynamical models to find the theoretical conditions under which measurements or samplings should be done to tell a null hypothesis $$H^0$$ from an alternative one $$H^a$$. Using the results of this model, we look at the time snapshots that increase the likelihood of obtaining significantly different results. We then test the model versus the experimental results.

As a proof of concept, we are testing Hyperion on an actual set of experiments. These experiments aim at understanding the environmental factors that could trigger phage induction on lysogens. The experiments are done with a *Salmonella* strain and four different species of phages. There are two variables that are being tested on the experiments:

1. Growth medium - Bacteria are cultured in two different media that control their growth rate and carrying capacity. In particular, bacteria cultured in Luria-Bertani medium (LB, for short) will have a higher growth rate and carrying capacity than bacteria cultured in Minimal Media (MM)

2. Bacterial strain - There are two bacterial strains: The wild type (WT) bacteria and the  *spoT-* (or mutants). It is expected that *spoT-* mutants are less inducible than WT.

The experiments are performed under the alternative hypothesis that the *spoT-* mutant is less likely to be induced than the WT. The null hypothesis is that the WT
grows slower than the *spoT-* mutant.

We show the results of the experiments performed. The tables show the results of the culturing in LB and MM media respectively. Notice that we only have final concentrations for the phages.


LB Media

||Bacterial concentration (cells/ml)|
|----|----------------|-------------------|---------------------|------------------------|                                                     
|Time (h)|WT |WT+CX |*spoT-* |*spoT-*+CX |
|0|$$2.20e7$$|$$2.20e7$$|$$1.80e7$$|$$1.80e7$$|
|2|$$2.93e7$$|$$2.93e7$$|$$3.30e7$$|$$3.30e7$$|
|5.5|$$2.93e8$$|$$2.93e8$$|$$1.47e8$$|$$1.47e8$$|
|11.5|$$4.50e8$$|$$4.00e8$$|$$3.33e8$$|$$1.93e8$$| 
|17|$$7.07e8$$|$$6.20e8$$|$$7.73e8$$|$$6.53e8$$|                              
|21|$$8.10e8$$|$$8.40e8$$|$$8.60e8$$|$$9.20e8$$|
|30|$$1.39e9$$|$$1.24e9$$|$$1.54e9$$|$$1.44e9$$| 


||Phage concentration (phages/ml)|
|----|----------------|-------------------|---------------------|------------------------|                   
|Time (h)|WT |WT+CX |*spoT-* |*spoT-*+CX |
|30|$$1.7e2$$|$$2.2e4$$|$$4e1$$|$$8e3$$|


MM

||Bacterial concentration (cells/ml)|
|----|----------------|-------------------|---------------------|------------------------|
|Time (h) |WT|WT+CX|*spoT-*|*spoT-*+CX|
|0|$$1.00e7$$|$$1.00e7$$|$$8.00e6$$|$$8.00e6$$|
|4|$$4.30e7$$|$$4.30e7$$|$$2.90e7$$|$$2.90e7$$|
|10|$$7.50e7$$|$$5.00e7$$|$$4.20e7$$|$$4.80e7$$|
|14|$$8.10e7$$|$$6.50e7$$|$$9.90e7$$|$$6.70e7$$|
|17|$$8.10e7$$|$$6.50e7$$|$$9.90e7$$|$$6.70e7$$|
|24|$$8.10e7$$|$$6.50e7$$|$$9.90e7$$|$$6.70e7$$|


||Phage concentration (phages/ml)|
|----|----------------|-------------------|---------------------|------------------------|                   
|Time (h)|WT |WT+CX |*spoT-* |*spoT-*+CX |
|24|$$1.7e2$$|$$2.2e4$$|$$4e1$$|$$8e3$$|