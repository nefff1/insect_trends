This repository contains code that was used in the analyses for the following manuscript:

Neff F, Korner-Nievergelt F, Rey E, Albrecht M, Bollmann K, Cahenzli F, Chittaro YD, Gossner MM, Martínez-Núñez C, Meier ES, Monnerat C, Moretti M, Roth T, Herzog F, Knop E. **Insect declines and increases linked more to climate than land-use changes**.

The following R code files are included in the folder *R_Code*:

-   **R_Main_analyses.R**: Main file containing all the relevant analyses reported in the study. Relies on processed data originating from analyses based on code in other files.

-   **R_Climate_change.R**: Code used to compute climate change variables

-   **R_Land_use_change.R**: Code used to compute land-use change variables

-   **R_Temperature_niche.R**: Code used to compute species' temperature niches

-   **R_Proportion_Switzerland.R**: Code used to calculate species distribution proportions in Switzerland

-   **R_Occupancy_detection_models.R**: Code used to prepare species records data and run occupancy-detection models. This code was run on a HPC cluster.

-   **f_occ_det.R**: Helper function used to run occupancy-detection models in Stan.

The following Stan code files are included in the folder *Stan_Code*:

-   **Stan_regression.stan**: Stan model code of the regression model fitted to species trends.

-   **Stan_occ_det_cmdstan.stan**: Stan model code of the occupancy-detection models.

The folder *Other* contains:

-   **Bash_start_HPC_cluster.sh**: example of the shell code used to run occupancy-detection models on the HPC cluster

-   **specieslist.txt**: List of all 390 study species and the insect groups they belong to
