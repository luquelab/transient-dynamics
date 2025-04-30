# Tutorial

## Description

This tutorial guides you through running a Jupyter notebook to simulate Lotka-Volterra predator-prey models (standard and adaptive). Clone this repository and run the notebook in Visual Studio Code to store generated data files. You can also use Google Colab, but you'll need to manually compare files with those in the GitHub repository.

## Running the Lotka-Volterra Jupyter notebook 

1. Clone this repository to your computer. This can be done via serveral methods: utilize GitHub Desktop app to clone the repository directly to your desktop, or `git clone` command on Anacoda prompt. 
<br>
<bold>Note</bold>: One can also open the notebook on the [Colab notebook](https://colab.research.google.com/github/luquelab/transient-dynamics/blob/main/notebooks/transient_dynamics.ipynb), but it is important to modify folder pathing in further steps.

2. Open the Jupyter notebook file (**transient_dynamics.ipynb**), follow sequentially and run the blocks sequentially until reaching the Directory Structure Setup, if the repository is cloned continue but if operating in Colab notebook change the parent directory to a folder name of choice for ease of finding. 

3. Continue to the model which you would like to simulate; the notebook contains a "Standard Lotka-Volterra Simulation" section and a "Lotka-Volterra with Carrying Capacity." Input parameters for the given variables as seen in the "Input parameter" section into the model parameter inputs section. If you would like to replicate the scenarios from the manuscript, open pathway (**theta_dynamics_code/data/transient-dynamics_output-current/example_scenarios**) and select a model (either standard_model or carrying_capacity_model) and then select a scenario. After which go to the **input.txt** file to copy parameter into the model parameter inputs section for the model selected then run this block. 

4. Run all blocks in the Analysis and Visualization section to produce figures and an Excel sheet labeled, **array_outputs.xslx**, with population, regime, weight, and processes data. Non-array variables, such as error analysis per regime, tipping points and initial parameters, will be stored in a text file labeled **non_array_outputs.txt**.

5. These files are stored in the **generated-data** folder (located by default in **theta_dynamics_code/data/transient-dynamics_output-current/data_output** if another folder or directory was not generated), containing the folders for both of the models: standard_model and carrying_capacity_model. Within these folder are subfolders for the figures and results labeled  **figures_standard** or **figures_carrying_capacity** and **results_standard** or **results_carrying_capacity**. The figures folder contains all graphs and plots made, and the results contains the excel and text files. 

6. If one utilized an **input.txt** file for a given scenario they can now compare the generated-data folder from simulation to the example_outputs folder scenarios.

## Additional notes

+ The same process described above can be used to generate simulation for both standard and carrying capacity models, one just needs to input parameters for both models and continue sequentially the upstream blocks. 

+ If you want to compare a scenario within another scenario in phase diagram, go to the "Scenario comparison simulation and analysis" block and input parameters in the "Input parameters of auxiliary scenario" block. Run all blocks to the Visualizing both scenarios within a single-phase diagram block. The figure will be saved into the figures_standard subfolder within the larger **generated-data** folder. 

+ Issues can occur if you try to run a block multiple times so ensure prior to usage that the outputs have been cleared before running the Notebook and following the steps. 