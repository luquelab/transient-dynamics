# General information
This file provides complementary information for the data file `Input_File_Kick_Test.csv`:

## Description of columns

+ `description`: short description of the experiment indicating active mechanisms, percentage by which concentrations of bacteria and/or phage have been kicked, and whether or not a logistic growth was implemented.

**Model parameters**
+ `r`: Bacterial growth rate constant. Units: 1/h.
+ `a`: Phage adsorption rate constant (known in code as 'd'). Units: ml/h.
+ `c`: Burst size. Units: Dimensionless.
+ `m`: Phage decay rate. Units: 1/h.
+ `K`: Bacteria carrying capacity. Units cells/ml.
+ `tau`: Observational time scale. Units: h.
+ `V`: Volume of the system. Units: ml.

**Predicted dynamic regime shift values**
+ `w_r`: Weight for the bacterial growth process. Units: Dimensionless.
+ `w_m`: Weight for the phage decay process. Units: Dimensionless.
+ `tau_r`: Critical timescale of the bacterial growth process. Units: h.
+ `tau_m`: Critical timescale of the phage decay process. Units: h.
+ `B_c`: Critical concentration of bacteria. Units: cells/ml.
+ `P_c`: Critical concentration of phages. Units: cells/ml.

**Initial conditions**
+ `B_0`: Initial concentration of bacteria. Units: cells/ml.
+ `P_0`: Initial concentration of phages. Units: virions/ml.
+ `t_f`: Final simulation time. Units: h.

**Experiment values**

+ `kick_b`: Percentage in which the concentration of bacteria is modified. Units: dimensionless.
+ `kick_p`: Percentage in which the concentration of phages is modified. Units: dimensionless.
+ `B0_kick`: Bacterial concentration after the kick. Units: cells/ml.
+ `P0_kick`: Phage concentration afther the kick. Units: virions/ml.
+ `t_kick`: Time at which the kick occurs. Units: h.

**Simulation values**
+ `B_f`: Final concentration of bacteria. Units: cells/ml.
+ `P_f`: Final concentration of phages. Units: virions/ml.
+ `a_B`: Amplitude of oscillations of bacteria (difference between maximum and minimum concentrations). Units: cells/ml.
+ `a_P`: Amplitude of oscillations of phage (difference between maximum and minimum concentrations). Units: virions/ml.

**Stability theory**
+ `B_eq`: Equilibrium concentration of bacteria. Units. cells/ml.
+ `P_eq`: Equilibrium concentration of phage. Units: virions/ml.

## Notes

**Initial concentrations**
+ Initial concentrations are always equilibrium concentrations in this case: `B_0`=`B_eq` and `P_0`=`P_eq`

**Implementation of the kick experiment**
+ The kick experiment can be performed in two ways: one is by a relative change of the initial concentrations, represented by `kick_b` and `kick_p`. Another one is to set absolute values `B0_kick` and `P0_kick`.

**Details of the numerical simulation**
+ The model was solved using the python library 'scipy' and the function 'solve_ivp'. The code used ran the Python3 version 3.6.9



