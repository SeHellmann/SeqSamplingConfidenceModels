---
---
---

# Modelling Reaction Time and Confidence Distributions in Decision Making

This repository contains data as well as code used in the paper **Simultaneous modeling of choice, confidence and response time in visual perception** (Hellmann, S., Zehetleitner, M., & Rausch, M. (in press, *Psychological Review*). 
[Preregistration on OSF](https://osf.io/x548k/) was developed and published with another data set from a previous study (Experiment 2 in Rausch, Hellmann, & Zehetleitner (2018). The other experiments are available in this repository. The most recent author manuscript of the article is available here: [https://osf.io/mrh78](https://osf.io/mrh78). 

## Structure:

-   dynWEV-source package file (.tar.gz)
-   dynConfiR-source package file (.tar.gz) (for the additional analyses for review). This is the 0.0.1 version of the package with likelihood and fitting functions, which is now available on [CRAN](https://cran.r-project.org/web/packages/dynConfiR/index.html). 
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
- The folder *Additional_Analyses* with further non-preregistered analyses conducted for the review process
    -   a script R-file *additional_analyses_for_review.R* with code to fit and predict two models (DDMConf and dynVis) and saved model fits for both models and both experiments
    -   two .RData files, *collected_fitsNpredicts_*Experiment*_review.RData* with all model fits and predictions for visualization
    -   a script R-file and folder for a small *parameter recovery* study for the dynWEV model
    -   *gen_model_weights.R* with code to transfer BIC values to model weights
    -   *simulation_dynWEV.R* with code to produce a figure with simulations for the dynWEV model with differen weight parameters
    -   *AUC_Tau_plot.R* with code to produce Supplementary Figure 1 for the relationship between metacognitive sensitivity and postdecisional accumulation time



## Usage:

-   Start R with package file in working directory

<!-- -->

    install.packages("dynWEV_0.0.tar.gz", type = "source", dependencies=TRUE,repos="http://a.cran.mirror")
    install.packages("dynWEV_0.0.1.tar.gz", type = "source", dependencies=TRUE,repos="http://a.cran.mirror")

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
