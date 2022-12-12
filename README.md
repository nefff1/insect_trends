Insect occupancy trends 1980-2020 in Switzerland
======

[![DOI](https://zenodo.org/badge/485801411.svg)](https://zenodo.org/badge/latestdoi/485801411)

This repository contains code that was used in the analyses for the following publication:

Neff F, Korner-Nievergelt F, Rey E, Albrecht M, Bollmann K, Cahenzli F, Chittaro Y, Gossner MM, Martínez-Núñez C, Meier ES, Monnerat C, Moretti M, Roth T, Herzog F, Knop E. 2022. **Different roles of concurring climate and regional land-use changes in past 40 years' insect trends**. Nature Communications. DOI: [10.1038/s41467-022-35223-3](https://doi.org/10.1038/s41467-022-35223-3).

The following R code files are included in the folder *R_Code*:

-   **R_Main_analyses.R**: Main file containing all the relevant analyses reported in the study. Relies on processed data originating from analyses based on code in other files.

-   **R_Climate_change.R**: Code used to compute climate change variables.

-   **R_Land_use_change.R**: Code used to compute land-use change variables.

-   **R_Temperature_niche.R**: Code used to compute species' temperature niches.

-   **R_Proportion_Switzerland.R**: Code used to calculate species distribution proportions in Switzerland.

-   **R_Occupancy_detection_models.R**: Code used to prepare species records data and run occupancy-detection models. This code was run on a HPC cluster.

-   **f_occ_det.R**: Helper function used to run occupancy-detection models in Stan.

The following Stan code files are included in the folder *Stan_Code*:

-   **Stan_regression_full.stan**: Stan model code of the regression model fitted to species trends. Full model used for versions 2 and 3.

-   **Stan_regression_restricted.stan**: Stan model code of the regression model fitted to species trends. Restricted model used for version 1 (parameter estimates for agricultural area and grassland-use intensity change only rely on species of agriculturally influenced habitats).

-   **Stan_occ_det_cmdstan.stan**: Stan model code of the occupancy-detection models.

The folder *Other* contains:

-   **Bash_start_HPC_cluster.sh**: example of the shell code used to run occupancy-detection models on the HPC cluster.

-   **specieslist.txt**: List of all 390 study species and the insect groups they belong to.

-   **Agricultural_species.txt**: List of species (partly) bound to agriculturally influenced habitats.
