---
---
---

# Modelling Reaction Time and Confidence Distributions in Decision Making

This repository contains data as well as code used in the paper **Simultaneous modeling of choice, confidence and response time in visual perception**. [Preregistration on OSF](https://osf.io/x548k/) was developed and published with another data set from a previous study (Experiment 2 in Rausch, Hellmann, & Zehetleitner (2018). The other experiments are available in this repository.

## Structure:

-   dynWEV-source package file (.tar.gz)
-   folders for the experiments analyzed in the study. The two motion discrimination datasets from experiment 2 are again included in sub-folders (and have individual files for the first three bullets). Following files are in the experiment folders:
    -   a *experiment* folder containing the files for running the experiments in Psychopy
    -   a .csv file ('data*Experiment*.csv') containing the raw data
    -   a script R-file for the actual analyses ('Script_FitNPredict_SeqSampConfModels\_*Experiment*.R'), including:
        -   Reading, preprocessing, and Aggregating Data
        -   Fitting model parameters
        -   Prediction of confidence and RT distributions and aggregation of predictions
    -   files to generate reported results, figures and tables in the paper (gen_descr_plots.R, gen_model_plots_and_BICAnalysis.R, and gen_table_fittedparameters.R)
    -   a script R-file for model identification analysis ('Script_ModelMimikryAnalysis.R')
    -   a *autosave_mimikry* folder with saved results from the model identification analysis
    -   a *saved_fits* folder with two files containing the fitted parameters from the experiment for diffusion based models ('fits_2DSD_WEV.RData') and race models ('fits_RacingModels.R'), respectively

## Usage:

-   Start R with package file in working directory

<!-- -->

    install.packages("dynWEV_0.0.tar.gz", type = "source", dependencies=TRUE,repos="http://a.cran.mirror")

-   If necessary, install required packages:

<!-- -->

    install.packages(c("plyr", "snow", "doSNOW", "BayesFactor", "tidyverse", "RColorBrewer", "gridExtra"))

-   To redo the whole analyses, run 'Script_FitNPredict_SeqSampConfModels\_*Experiment*.R' in the respective experiment folder
-   Use the saved results in 'collected_fitsNpredicts.RData' in the respective experiment folders for all other analyses

## Compatibility for package versions

As some R packages are under constant development we included the file sessionInfo.txt with the necessary information about versions of R packages used for the original analyses.

### References

Rausch, M., Hellmann, S. & Zehetleitner, M. Confidence in masked orientation judgments is informed by both evidence and visibility. Atten Percept Psychophys 80, 134-154 (2018). <https://doi.org/10.3758/s13414-017-1431-5>

## Contact

For comments, remarks, and questions please contact me: [sebastian.hellmann\@ku.de](mailto:sebastian.hellmann@ku.de){.email}
