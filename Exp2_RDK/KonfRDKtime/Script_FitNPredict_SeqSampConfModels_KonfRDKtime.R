###########################################################################
#####      Comparison Sequential Sampling Confidence Models        ########
###########################################################################

# Sebastian Hellmann, 07.05.2021

# 1) Read in data, preprocess and aggregate data for later visualization
# 2) Fit the models (dynWEV, 2DSD, IRM(t) and PCRM(t))  
# 3) predict rating and rt distribution and aggregate for visualization


###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/AccumulatorModels")
# or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
{
  library(plyr)
  library(snow)
  library(doSNOW)
  library(tidyverse)
  library(dynWEV)
}
#=============================================================
#=============================================================
###   Load results from previous analysis   ####
### (Includes results from sections 1) & 2) ###
if (file.exists("collected_fitsNpredicts.RData")) {
  load("collected_fitsNpredicts.RData")
} else {
#=============================================================
#=============================================================
########### 1) Read, Preprocess, Aggregate Data ##############
####  Read in data and exclude trials and participants   #####
Data <- read.table("dataKonfRDKTime.csv", header=TRUE, sep=",")
Data <- group_by(Data, participant)

Ntotal <-  Data %>%
  summarise(Ntot = sum(n())) 
#### Sample descriptives         ####
print("Gender distribution:")
print(Data %>% distinct(participant, gender) %>% group_by(gender) %>% summarise(N = n()))
print("Age distribution:")
print(Data %>% distinct(participant, age) %>% 
        group_by(participant) %>%
        summarise(age=mean(age)) %>%  #Because one participant got older over the sessions, we take the mean age of participants
        ungroup() %>% 
        summarise(min = min(age), max= max(age), M = mean(age), SD=sd(age)))


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

CheckParticipants <- Data %>%
  summarise(
    nTrials = length(correct),
    Performance = mean(correct),
    AboveChance = as.numeric(try(extractBF(proportionBF(y=sum(correct), N=length(correct), p=1/2))$bf)),
    WhichModeRating = Mode(Rating),
    NumModeRating = sum(Rating  == Mode(Rating))/length(Rating),
    nSessions = length(unique(session)))

BadSubjects <- CheckParticipants %>%
  filter(Performance < .5 | AboveChance < 3 | 
           NumModeRating > .9) # | nSessions < 3

BadSubjects

Data <- subset(Data, !participant %in% BadSubjects$participant)

Data <- Data %>%
  filter(rt < mean(rt) + 4*sd(rt)  & rt > .3)


Data <- Data %>%
  filter(min(table(coherence, stimulus)) >= 20 )

Nanalysis <-  left_join(summarise(Data, Nana = sum(n())), Ntotal)
cat(paste("Excluded Participants because of missing sessions:", paste((subset(BadSubjects, nSessions < 3))$participant, collapse = ", "), "\n"))
cat(paste("Excluded Participants based on Accuracy:", paste((subset(BadSubjects, nSessions == 3))$participant, collapse = ", "), "\n"))
cat("Summary of Proportion of Excluded Trials per (non-excluded) Participant :\n")
print(summary(1- Nanalysis$Nana/Nanalysis$Ntot))

Data <- rename(Data, condition=coherence)
nRatings <- length(unique(Data$rating))
nConds <- length(unique(Data$condition))
cond_levels <- sort(unique(Data$condition))
Data <- Data %>% 
  mutate(condition = factor(condition,levels = cond_levels, 
                                           labels = as.character(cond_levels*100)))

#### Aggregate Data                                        #####
#### Compute confidence rating distribution of the Data    #####
Data <- Data %>% group_by(stimulus, condition, participant) %>%
  mutate(nrows=n())
Data_RatingDist_part_stim <- Data %>% 
  group_by(stimulus, condition, rating, response, correct, participant) %>% 
  summarise(p = n()/(mean(nrows))) %>% ungroup() %>%
  full_join(y = expand.grid(stimulus = unique(Data$stimulus), 
                          response = unique(Data$stimulus), 
                          condition = unique(Data$condition),
                          rating = 1:nRatings, 
                          participant = unique(Data$participant))) %>%
  mutate(p = ifelse(is.na(p), 0, p), correct = as.numeric(response == stimulus)) 
Data_RatingDist_part <- Data_RatingDist_part_stim %>%
  group_by(correct, condition, rating, participant) %>% 
  summarise(p = mean(p))
### Sanity Checks:
# sum(Data_RatingDist_part$p)
# table((Data_RatingDist_part %>% group_by(condition, participant) %>% summarise(p = sum(p)))$p)

# For the plots we won't differentiate between 
# stimulus directions and participants
Data_RatingDist_corr_cond <- Data_RatingDist_part %>% 
  group_by(correct, condition, rating) %>% 
  summarise(p = mean(p))
#sum(Data_RatingDist_corr_cond$p)

#### Compute Mean Rating of the Data    #####
Data_MRating_corr_cond <- Data_RatingDist_part %>% 
  group_by(condition, participant, correct) %>%
  summarise(MRating = sum(rating*p)/sum(p)) %>%
  group_by(condition, correct) %>%
  summarise(SER = sd(MRating, na.rm = TRUE)/sqrt(n()),
            MRating = mean(MRating, na.rm = TRUE))

#### Compute Reaction Time Quantiles of the Data grouped by rating and accuracy  #####
Data_RTQuants_corr_rating <- Data %>%
  group_by(rating, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 


####################### 2) Fit Models and Predict Distributions   #################################
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores, "SOCK", outfile = "")
registerDoSNOW(cl)
clusterEvalQ(cl, library(dynWEV))

######### Start fitting with models including postdecisional accumulation of evidence #############
fitData <- cbind(rbind(Data, Data), model = rep(c("2DSD", "WEVmu"), each = nrow(Data)))
clusterExport(cl, "fitData")

t00 <- Sys.time()
fits_WEVmodels <- ddply(fitData,.(participant, model),
                              function(df) fitWEV(df, df$model[1], logging = TRUE, 
                                                  restr_tau ="simult_conf", nRatings = 5),
                              .parallel = TRUE) 
# we parallelize over participants, here, 
# and not within the fitting process

dir.create("saved_fits")
save(file = "saved_fits/fits_2DSD_WEV.RData", fits_WEVmodels)
print(paste("Fitting 2DSD and WEV took...",
            as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
            " mins"))

##### Fit Race Models/ bounded accumulation models  #######################
fitData <- cbind(rbind(Data, Data, Data, Data), 
                 model = rep(c("IRM", "PCRM", "IRMt", "PCRMt"), each = nrow(Data)))
clusterExport(cl, "fitData")

t00 <- Sys.time()
fits_RMmodels <- ddply(fitData,.(participant, model),
                       function(df) fitRM(df, df$model[1], time_scaled = FALSE,
                                          logging = TRUE, nRatings = 5),
                       .parallel = TRUE)
dir.create("saved_fits")
save(file = "saved_fits/fits_RacingModels.RData", fits_RMmodels)
print(paste("Fitting both Racing Models took...",
            as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
            " mins"))

#====================================================================================
########### 3) Compute model predictions with fitted parameters #####################
#####   Prediction of response and rating distributions  #####
clusterExport(cl, c("fits_WEVmodels"))
clusterExport(cl, c("fits_RMmodels"))
maxrt <- max(Data$rt, 20)
clusterExport(cl, c("maxrt"))
preds_WEV <- ddply(fits_WEVmodels,.(participant, model),
                   function(df) predictWEV_Conf(df, model = df$model[1], subdivisions = 1000, 
                                                simult_conf=TRUE, maxrt = maxrt),
                   .parallel = T)
table(preds_WEV$info)
preds_RM <- ddply(fits_RMmodels,.(participant, model),
                  function(df) predictRM_Conf(df, model = df$model[1], subdivisions = 1000, time_scaled = FALSE, maxrt = maxrt),
                  .parallel = TRUE)
table(preds_RM$info)


#####          Prediction of RT densities              #######
maxrt <- max(Data$rt, 6)   ## For visualization a maximum RT of 6 is enough
compute_RTdens_from_fits <- function(Confpred_data) {
  if (Confpred_data$model[1] %in% c("IRM","PCRM","IRMt", "PCRMt")) {
    paramDf <- subset(fits_RMmodels, 
                      model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
    res <- predictRM_RT(paramDf = paramDf, 
                        model = Confpred_data$model[1], 
                        maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0),
                        scaled = TRUE, DistConf = Confpred_data,
                        .progress = FALSE)
  } else {
    paramDf <- subset(fits_WEVmodels, 
                      model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
    res <- predictWEV_RT(paramDf = paramDf, 
                         model = Confpred_data$model[1], 
                         maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0), 
                         scaled = TRUE, simult_conf=TRUE, DistConf = Confpred_data,
                         .progress = FALSE)
  }
  return(res)
}

clusterExport(cl, c("preds_WEV", "preds_RM", "maxrt", "compute_RTdens_from_fits"))
RT_dist_RM <- ddply(preds_RM,.(participant, model),
                    .fun = compute_RTdens_from_fits, 
                    .parallel = T)
RT_dist_RM <- RT_dist_RM %>% group_by(model, participant, correct, rating, condition, rt) %>%
  summarise(dens = mean(dens), densscaled = mean(densscaled))

RT_dist_WEV <- ddply(preds_WEV,.(participant, model),
                     .fun = compute_RTdens_from_fits, 
                     .parallel = T)
RT_dist_WEV <- RT_dist_WEV %>% group_by(model, participant, correct, rating, condition, rt) %>%
  summarise(dens = mean(dens), densscaled = mean(densscaled))
RT_dist <- rbind(RT_dist_RM, RT_dist_WEV)

stopCluster(cl)
rm(list = c("cl", "RT_dist_WEV", "RT_dist_RM"))


#=============================================================
#####   Aggregating Predictions for Visualization     ########
# For the visualization we aggregate always over stimulus and
# response category, because we are only interested in correct
# vs. incorrect responses;
# Aggregation is always within each model (keeping them distinct)

##   Combine and Aggregate Confidence Rating Distribution   ##
Preds_RatingDist_corr_cond <- rbind(preds_WEV, preds_RM) %>%
  group_by(model, rating, correct, condition) %>%
  summarise(p = mean(p)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")), 
         condition = factor(condition,levels = 1:nConds, 
                            labels = as.character(cond_levels*100)))


# # # Sanity checks:
# Preds_RatingDist_corr_cond %>% filter(model %in% c("IRM", "IRMt", "PCRM", "PCRMt")) %>%
#   group_by(model, condition) %>% summarise(p = sum(p))
# Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
#   summarise(p=sum(p))


##    Compute Mean Rating Accross Conditions 
# This is good to visualize a folded-X- and double-increase-pattern 
Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% group_by(model, condition, correct) %>%
  summarise(MRating = sum(p*rating)/sum(p))

##   Combine and Aggregate RT densities 
Ns_part <- Data %>% group_by(participant) %>% summarise(N=n(), MinRT = min(rt))  %>%
  select(participant, N)
Preds_RTdens_corr_cond_rating <- RT_dist %>% 
  left_join(Ns_part) %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data)) %>%  
  # Use a weighted mean, here (only affects visualisation), because participants had different number of trials 
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")), 
         condition = factor(condition,levels = 1:nConds, 
                            labels = as.character(cond_levels*100)))  


####    Computation of RT quantiles    
Preds_RTQuants_corr_rating <- Preds_RTdens_corr_cond_rating %>% group_by(model, rt, correct, rating) %>%
  summarise(dens = mean(dens))%>%
  mutate(model = factor(model, levels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels=c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"))) %>%
  group_by(model, correct, rating) %>% 
  do(RTDensityToQuantiles(.)) %>% 
  left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                               model, rating, correct), 
                      p_correct=mean(p)))

#=============================================================
#####          Clean Up and Save Results           ###########
gc()

save(file="collected_fitsNpredicts.RData", 
     Data , nConds, nRatings,  maxrt, cond_levels,
     Data_RatingDist_part ,  Data_RatingDist_corr_cond ,  
     Data_MRating_corr_cond ,
     fits_WEVmodels ,  fits_RMmodels ,
     preds_RM ,  preds_WEV , Ns_part,
     Preds_RatingDist_corr_cond ,  Preds_MRating_corr_cond , 
     Preds_RTdens_corr_cond_rating  , 
     Preds_RTQuants_corr_rating, Data_RTQuants_corr_rating)

print("End of fitting and prediction. Saved all results.")
}
###############################################################

###### REST ##################
#rm(list=ls())
