# Luque Lab project template

## WHAT
This is a tool that follows the protocol for FODAM (Finite Observational Dynamic Analysis Method) which was first proposed in the Manuscript Emerging dynamic regimes and tipping points from finite empirical principles by Sergio Cobo-Lopez, Matthew Witt, Forest L. Rohwer, and Antoni Luque (2023). This tool is made independently from the manuscript's tool, using a conditional variable (theta) that is either equivalent to one under the condition that it's corresponding process weight is above the threshold (1) at a given time or zero if the weight is below the threshold. This allows the tool to eliminate non-significant processes at given points producing a adaptive model comparable to that of the full model (in this case the Lotka-Volterra Predator-Prey model).

With this tool we hope to further validate the findings of the manuscript and allow for a more comprehensive interpretation of dynamical models, one that would allow one to make both predictions and understand fully what variables are most pertinent in causing shifts in agents levels. We will continue to update this tool and expand it to other dynamic systems such as the SIR model, Lotka-Volterra Competitor model, and exploring phage decay.

## WHO
The template was originally concieved by Antoni Luque based on recommendations from Noble PLoS Comp Biol 2009, Wilson et al PLoS Comp Biol 2017, and Briney "Data Management for Researchers" (2015), Hunt and Thomas "The Pragmatic Programmer: Your Journey to Mastery" (2019 2nd ed).

## CONTRIBUTORS
This tool was produced by @luquelab, inspired by the original manuscript's repository (needle-finder) by Sergio Cobo-López. The main contributors are Lucas J. Carbajal and Antoni Luque, with additional assistance and guidance provided by the rest of the members in the Luque Lab.

## WHEN
This is an evolving repository Started: 2024-09-11

End: Ongoing

## FILES & FOLDERS
FOLDER: /bin
--> This folder contains basic scripts and executable files.

FOLDER: /data
--> This folder contains the raw data associated with the project and the potential references.

FOLDER: /doc
--> This folder contains the manuscript, digital copies of the cited references, figures, and other associated files for publication.

FOLDER: /results
--> This folder contains the results, performance analysis, and commented references associated with the project.

FOLDER: /src
--> This folder contains source code files that require compilation.

The syntax of markdown files (.md) is CommonMark unless specified otherwise (https://commonmark.org/help/)
