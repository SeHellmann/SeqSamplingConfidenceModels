---
---
---

# Modelling Reaction Time and Confidence Distributions in Decision Making

This repository contains data as well as code used in the paper **Simultaneous modeling of choice, confidence and response time in visual perception** (Hellmann, S., Zehetleitner, M., & Rausch, M., 2023, *Psychological Review*, [doi: 10.1037/rev0000411](https://doi.org/10.1037/rev0000411)). 
[Preregistration on OSF](https://osf.io/x548k/) was developed and published with another data set from a previous study (Experiment 2 in Rausch, Hellmann, & Zehetleitner (2018). The other experiments are available in this repository. The most recent author manuscript of the article is available here: [https://osf.io/mzfkr](https://osf.io/mzfkr). 

## Structure:

-   dynWEV-source package file (.tar.gz)
-   dynConfiR-source package file (.tar.gz) (for the additional analyses for review). This is the 0.0.1 version of the package with likelihood and fitting functions, which is now available on [GitHub](https://github.com/SeHellmann/dynConfiR). 
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

### Setup the environment:

-   Start with a completely new R 4.0.5 installation
-   On windows, install rtools40 (https://cran.r-project.org/bin/windows/Rtools/rtools40.html)
-   Install the renv package using `install.packages("renv")`
-   Change the working directory to the project directory `setwd('~/Material_for_SeqSamplingModelsOfChoiceConfRT')`
-   Use `renv::restore('renv.lock')` to install all packages with their respective version. 
    Note, that this will install the packages in your default library of your R-4.0.5 installation!
-   Install the local packages:
```
    install.packages('dynWEV_0.0.tar.gz', repos = NULL, type = 'source')
    install.packages('dynConfiR_0.0.1.tar.gz', repos=NULL, type = 'source')
```
### Redo whole analysis:
-   To redo the whole analyses with modelling fitting:
    -   remove the files with the results 'collected_fitsNpredicts.RData' in the respective experiment folder
    -   run 'Script_FitNPredict_SeqSampConfModels\_*Experiment*.R' from within the respective experiment folder
-   To reuse the quantitative model comparison and produce the figures, using the already computed model fits, simply use the saved results in 'collected_fitsNpredicts.RData' in the respective experiment folders for all other analyses
-   Note that the scripts 'Script_ModelMimikryAnalysis.R' always run the recovery analysis, irrespective of whether saved results are present or not and that this may take considerable time!

<!--

## Compatibility for package versions

As some R packages are under constant development we included the file sessionInfo.txt with the necessary information about versions of R packages used for the original analyses.
 -->
 
### References

Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence, and response time in visual perception. Psychological Review. Advance online publication. <https://doi.org/10.1037/rev0000411>


## Contact

For comments, remarks, and questions please contact me: [sebastian.hellmann\@tum.de](mailto:sebastian.hellmann@tum.de){.email}
