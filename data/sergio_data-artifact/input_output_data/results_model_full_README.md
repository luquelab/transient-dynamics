# General information
This file provides complementary information for the data file `values_r_active.csv`:

## Description of columns

**Model parameters**
+ `r`: Bacterial growth rate constant. Units: 1/h.
+ `a`: Phage adsorption rate constant. Units: ml/h.
+ `c`: Burst size. Units: Dimensionless.
+ `m`: Phage decay rate. Units: 1/h.
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

**Simulation values**
+ `B_f`: Final concentration of bacteria. Units: cells/ml.
+ `P_f`: Final concentration of phages. Units: virions/ml.
+ `t_B_V`: Time at which there is one bacterium in the system's volume. Units: h.
+ `t_P_V`: Time at which there is one phage in the system's volume. Units: h.
+ `t_B_c_1`: Time at which the bacterial concentration is equal to the critical concentration (first time). Units: h.
+ `t_P_c_1`: Time at which the bacterial concentration is equal to the critical concentration (first time). Units: h.
+ `t_B_c_2`: Time at which the bacterial concentration is equal to the critical concentration (second time). Units: h.
+ `t_P_c_2`: Time at which the bacterial concentration is equal to the critical concentration (second time). Units: h.

**Stability theory**
+ `B_eq`: Equilibrium concentration of bacteria. Units. cells/ml.
+ `P_eq`: Equilibrium concentration of phage. Units: virions/ml.
+ `lambda_1`: First eigenvalue for the asymptotic equilibrium. Units: 1/h.
+ `lambda_2`: Second eigenvalue for the asymptotic equilibrium. Units: 1/h.

## Notes
**Full model**
+ dB/dt = r*B - a*B*P
+ dP/dt = c*a*B*P - m*P

**Details of the numerical simulation**
+ The model was solved using the python library 'scipy' and the function 'solve_ivp'. The code used ran the Python3 version 3.6.9

**Stability theory analysis**
+ The stability was done for the full model.
+ The rates were set to zero to obtain equilibrium values.
+ The Jacobian was diagonalized to obtain the eigenvalues.
+ Stable solutions or centers were selected.
