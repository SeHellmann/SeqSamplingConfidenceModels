###########################################################################
#####      Additional analyses for the review process for           #######
#####           Sequential Sampling Confidence Models               #######
###########################################################################

# Sebastian Hellmann, 28.09.2022

# 1) Read in data, preprocess and aggregate data for later visualization
# 2) Fit the additional models (dynVib and DDMConf)
#     for both experiments
# 3) Predict rating and RT distribution and aggregate for visualization



###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/SeqSamplingConfidence")
# # or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
{
  library(Rmisc)
  library(plyr)
  library(snow)
  library(doSNOW)
  library(tidyverse)
  library(dynWEV)
  library(dynConfiR)
}
#=============================================================
#=============================================================
###   Fit the dynVib model to the KonfMaskTime data   ####

load("../Exp1_KonfMaskTime/collected_fitsNpredicts.RData")
print("Loaded data and previously computed fits and predictions from KonfMaskTime.")

fit_dynVib <- fitRTConfModels(Data, "dynWEV", nRatings = 5, fixed = list(sym_thetas=FALSE, w=0, tau=0),
                              restr_tau = "simult_conf", parallel="models", logging=TRUE)
pred_dynVib <- predictConfModels(fit_dynVib, simult_conf = TRUE, .progress = TRUE)
maxrt <- max(Data$rt, 6)
minrt <- min(fits_RMmodels$t0, fits_WEVmodels$t0)
predRT_dynVib <- predictRTModels(fit_dynVib, maxrt = maxrt, minrt = minrt,
                                 simult_conf = TRUE, scaled = TRUE,
                                 DistConf = pred_dynVib)
fit_dynVib <- rename(fit_dynVib, participant=sbj)
fit_dynVib$model <- "dynVib"
pred_dynVib <- rename(pred_dynVib, participant=sbj)
pred_dynVib$model <- "dynVib"
predRT_dynVib <- rename(predRT_dynVib, participant=sbj)
predRT_dynVib$model <- "dynVib"

save(fit_dynVib, pred_dynVib, predRT_dynVib,
     file="fitsNpredicts_dynVib_KonfMaskTime.RData")



###   Fit the DDMConf model to the KonfMaskTime data   ####

fit_DDMConf <- fitRTConfModels(Data, "DDMConf", nRatings = 5, fixed = list(sym_thetas=FALSE),
                               restr_tau = "simult_conf", parallel="models", logging=TRUE)
pred_DDMConf <- predictConfModels(fit_DDMConf, simult_conf = TRUE, .progress = TRUE)
maxrt <- max(Data$rt, 6)
minrt <- min(fits_RMmodels$t0, fits_WEVmodels$t0)
predRT_DDMConf <- predictRTModels(fit_DDMConf, maxrt = maxrt, minrt = minrt,
                                  simult_conf = TRUE, scaled = TRUE,
                                  DistConf = pred_DDMConf)
pred_DDMConf$correct <- as.numeric(if_else(pred_DDMConf$response=="upper", 1,-1)==pred_DDMConf$stimulus)
predRT_DDMConf$correct <- as.numeric(if_else(predRT_DDMConf$response=="upper", 1,-1)==predRT_DDMConf$stimulus)

save(fit_DDMConf, pred_DDMConf, predRT_DDMConf,
     file="fitsNpredicts_DDMConf_KonfMaskTime.RData")

#=============================================================
#####   Aggregating Predictions for Visualization     ########
load("../Exp1_KonfMaskTime/collected_fitsNpredicts.RData")
load("fitsNpredicts_DDMConf_KonfMaskTime.RData")
load("fitsNpredicts_dynVib_KonfMaskTime.RData")
# For the visualization we aggregate always over stimulus and
# response category, because we are only interested in correct
# vs. incorrect responses;
# Aggregation is always within each model (keeping them distinct)

##   Combine and Aggregate Confidence Rating Distribution   ##
Preds_RatingDist_corr_cond <- rbind(preds_WEV, preds_RM,
                                    pred_dynVib[,names(preds_RM)],
                                    pred_DDMConf[,names(preds_RM)]) %>%
  group_by(model, rating, correct, condition) %>%
  summarise(p = mean(p)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = paste0(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))
# # # Sanity checks:
# Preds_RatingDist_corr_cond %>% filter(model %in% c("IRM", "IRMt", "PCRM", "PCRMt")) %>%
#   group_by(model, condition) %>% summarise(p = sum(p))
# Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
#   summarise(p=sum(p))
Preds_RatingDist_corr_cond %>% filter(model=="DDMConf" & correct==1)

##    Compute Mean Rating Accross Conditions
# This is good to visualize a folded-X- and double-increase-pattern
Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% group_by(model, condition, correct) %>%
  summarise(MRating = sum(p*rating)/sum(p))
Preds_MRating_corr_cond <- rbind(preds_WEV, preds_RM,
                                 pred_dynVib[,names(preds_RM)],
                                 pred_DDMConf[,names(preds_RM)]) %>% 
  group_by(model, participant, condition, correct) %>%
  summarise(MRating = 0.5*sum(p*rating)/sum(p))%>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = paste0(c("8.3", "16.7", "33.3", "66.7", "133.3"), ""))) %>%
  group_by(model) %>%
  summarise(summarySEwithin(cur_data(),measurevar = "MRating", idvar = "participant", 
                  withinvars = c("condition", "correct"))) %>%
  rename(sewithin = se) %>% select(-MRating) %>%
right_join(mutate(Preds_MRating_corr_cond, correct=as.factor(correct)), by=c("condition", "correct", "model"))

##   Combine and Aggregate RT densities
Ns_part <- Data %>% group_by(participant) %>% summarise(N=n(), MinRT = min(rt))  %>%
  select(participant, N)

Preds_RTdens_corr_cond_rating_dynVib <- predRT_dynVib %>%
  left_join(Ns_part) %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = paste0(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))


Preds_RTdens_corr_cond_rating_DDMConf <- predRT_DDMConf %>%
  left_join(Ns_part) %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = paste0(c("8.3", "16.7", "33.3", "66.7", "133.3"), "")))

Preds_RTQuants_corr_rating <- rbind(Preds_RTdens_corr_cond_rating,
                                    Preds_RTdens_corr_cond_rating_dynVib,
                                    Preds_RTdens_corr_cond_rating_DDMConf) %>%
  group_by(model, rt, correct, rating) %>%
  summarise(dens = mean(dens)) %>%
  group_by(model, correct, rating) %>%
  do(RTDensityToQuantiles(.)) %>%
  left_join(summarise(group_by(Preds_RatingDist_corr_cond,
                               model, rating, correct),
                      p_correct=mean(p)))

Preds_RTQuants_corr_cond <- rbind(Preds_RTdens_corr_cond_rating,
                                    Preds_RTdens_corr_cond_rating_dynVib,
                                    Preds_RTdens_corr_cond_rating_DDMConf) %>%
  group_by(model, rt, correct, condition) %>%
  summarise(dens = mean(dens)) %>%
  group_by(model, correct, condition) %>%
  do(RTDensityToQuantiles(.)) %>%
  left_join(summarise(group_by(Preds_RatingDist_corr_cond,
                               model, condition, correct),
                      p_correct=sum(p)))

Data_RTQuants_corr_cond <- Data  %>%
  group_by(condition, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(Data_RatingDist_corr_cond, p_correct=sum(p)), by=c("correct","condition"))

#=============================================================
#####          Clean Up and Save Results           ###########
gc()

save(file="collected_fitsNpredicts_KonfMaskTime_review2.RData",
     Data , nConds, nRatings,  maxrt, cond_levels,
     Data_RatingDist_part ,  Data_RatingDist_corr_cond ,
     Data_MRating_corr_cond ,
     fits_WEVmodels ,  fits_RMmodels , fit_dynVib,fit_DDMConf,
     preds_RM ,  preds_WEV , pred_dynVib, pred_DDMConf, Ns_part,
     Preds_RatingDist_corr_cond ,  Preds_MRating_corr_cond ,
     Preds_RTdens_corr_cond_rating  ,
     Preds_RTQuants_corr_rating, Data_RTQuants_corr_rating,
     Preds_RTQuants_corr_cond , Data_RTQuants_corr_cond )

print("End of fitting and prediction. Saved all results.")

###############################################################












###   Fit the dynVib model to the RDK data   ####
load("../Exp2_RDK/collected_fitsNpredicts_bothRDKExperiments.RData")
print("Loaded data and previously computed fits and predictions from RDK experiments.")

fit_dynVib <- fitRTConfModels(Data, "dynWEV", nRatings = 5, fixed = list(sym_thetas=FALSE, w=0, tau=0),
                              restr_tau = "simult_conf", parallel="models", logging=TRUE)
pred_dynVib <- predictConfModels(fit_dynVib, simult_conf = TRUE, .progress = TRUE)
maxrt <- max(Data$rt, 6)
minrt <- min(fits_RMmodels$t0, fits_WEVmodels$t0)
predRT_dynVib <- predictRTModels(fit_dynVib, maxrt = maxrt, minrt = minrt,
                                 simult_conf = TRUE, scaled = TRUE,
                                 DistConf = pred_dynVib)
fit_dynVib$model <- "dynVib"
pred_dynVib$model <- "dynVib"
predRT_dynVib$model <- "dynVib"

save(fit_dynVib, pred_dynVib, predRT_dynVib,
     file="fitsNpredicts_dynVib_RDKexperiment.RData")

###   Fit the DDMConf model to the KonfMaskTime data   ####

fit_DDMConf <- fitRTConfModels(Data, "DDMConf", nRatings = 5, fixed = list(sym_thetas=FALSE),
                               restr_tau = "simult_conf", parallel="models", logging=TRUE)
pred_DDMConf <- predictConfModels(filter(fit_DDMConf,!is.na(negLogLik)), simult_conf = TRUE, .progress = TRUE)
maxrt <- max(Data$rt, 6)
minrt <- min(fits_RMmodels$t0, fits_WEVmodels$t0)
predRT_DDMConf <- predictRTModels(filter(fit_DDMConf,!is.na(negLogLik)), 
                                  maxrt = maxrt, minrt = minrt,
                                  simult_conf = TRUE, scaled = TRUE,
                                  DistConf = pred_DDMConf)
#pred_DDMConf$correct <- as.numeric(if_else(pred_DDMConf$response=="upper", 1,-1)==pred_DDMConf$stimulus)
#predRT_DDMConf$correct <- as.numeric(if_else(predRT_DDMConf$response=="upper", 1,-1)==predRT_DDMConf$stimulus)

save(fit_DDMConf, pred_DDMConf, predRT_DDMConf,
     file="fitsNpredicts_DDMConf_RDKexperiment.RData")



#=============================================================
#####   Aggregating Predictions for Visualization     ########
load("../Exp2_RDK/collected_fitsNpredicts_bothRDKExperiments.RData")
load("fitsNpredicts_DDMConf_RDKexperiment.RData")
load("fitsNpredicts_dynVib_RDKexperiment.RData")
# For the visualization we aggregate always over stimulus and
# response category, because we are only interested in correct
# vs. incorrect responses;
# Aggregation is always within each model (keeping them distinct)
cond_levels <- levels(Data$condition)

##   Combine and Aggregate Confidence Rating Distribution   ##
Preds_RatingDist_corr_cond <- rbind(preds_WEV, preds_RM,
                                    pred_dynVib[,names(preds_RM)],
                                    pred_DDMConf[,names(preds_RM)]) %>%
  group_by(model, rating, correct, condition) %>%
  summarise(p = mean(p)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = cond_levels))
# # # Sanity checks:
# Preds_RatingDist_corr_cond %>% filter(model %in% c("IRM", "IRMt", "PCRM", "PCRMt")) %>%
#   group_by(model, condition) %>% summarise(p = sum(p))
# Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
#   summarise(p=sum(p))
Preds_RatingDist_corr_cond %>% filter(model=="DDMConf" & correct==1)

##    Compute Mean Rating Accross Conditions
# This is good to visualize a folded-X- and double-increase-pattern
Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% group_by(model, condition, correct) %>%
  summarise(MRating = sum(p*rating)/sum(p))
Preds_MRating_corr_cond <- rbind(preds_WEV, preds_RM,
                                 pred_dynVib[,names(preds_RM)],
                                 pred_DDMConf[,names(preds_RM)]) %>% 
  group_by(model, participant, condition, correct) %>%
  summarise(MRating = 0.5*sum(p*rating)/sum(p))%>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = cond_levels)) %>%
  group_by(model) %>%
  summarise(summarySEwithin(cur_data(),measurevar = "MRating", idvar = "participant", 
                            withinvars = c("condition", "correct"))) %>%
  rename(sewithin = se) %>% select(-MRating) %>%
  right_join(mutate(Preds_MRating_corr_cond, correct=as.factor(correct)), by=c("condition", "correct", "model"))


##   Combine and Aggregate RT densities
Ns_part <- Data %>% group_by(participant) %>% summarise(N=n(), MinRT = min(rt))  %>%
  select(participant, N)

Preds_RTdens_corr_cond_rating_dynVib <- predRT_dynVib %>%
  left_join(Ns_part) %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = cond_levels))
nrow(Data)

# Preds_RTdens_corr_cond_rating_DDMConf <- predRT_DDMConf %>%
#   left_join(Ns_part) %>%
#   ungroup() %>%
#   select(-participant) %>%
#   group_by(rating, condition, model, correct, rt) %>%
#   summarise(nrows = n(),
#             dens = sum(dens*N)/nrow(subset(Data, participant %in% unique(pred_DDMConf$participant))), 
#             densscaled = sum(N*densscaled)/nrow(Data)) %>%
#   mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
#                                           "dynVib", "DDMConf"),
#                         labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
#                                    "dynVib", "DDMConf")),
#          condition = factor(condition,levels = 1:nConds,
#                             labels = cond_levels))
# 
Preds_RTdens_corr_cond_rating_DDMConf <- predRT_DDMConf %>%
  left_join(Ns_part) %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                          "dynVib", "DDMConf"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM",
                                   "dynVib", "DDMConf")),
         condition = factor(condition,levels = 1:nConds,
                            labels = cond_levels))

Preds_RTQuants_corr_rating <- rbind(Preds_RTdens_corr_cond_rating,
                                    Preds_RTdens_corr_cond_rating_dynVib,
                                    Preds_RTdens_corr_cond_rating_DDMConf) %>%
  group_by(model, rt, correct, rating) %>%
  summarise(dens = mean(dens)) %>%
  group_by(model, correct, rating) %>%
  do(RTDensityToQuantiles(.)) %>%
  left_join(summarise(group_by(Preds_RatingDist_corr_cond,
                               model, rating, correct),
                      p_correct=mean(p)))
Preds_RTQuants_corr_cond <- rbind(Preds_RTdens_corr_cond_rating,
                                  Preds_RTdens_corr_cond_rating_dynVib,
                                  Preds_RTdens_corr_cond_rating_DDMConf) %>%
  group_by(model, rt, correct, condition) %>%
  summarise(dens = mean(dens)) %>%
  group_by(model, correct, condition) %>%
  do(RTDensityToQuantiles(.)) %>%
  left_join(summarise(group_by(Preds_RatingDist_corr_cond,
                               model, condition, correct),
                      p_correct=sum(p)))

Data_RTQuants_corr_cond <- Data  %>%
  group_by(condition, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
  left_join(summarise(Data_RatingDist_corr_cond, p_correct=sum(p)), by=c("correct","condition"))

#=============================================================
#####          Clean Up and Save Results           ###########
gc()

save(file="collected_fitsNpredicts_RDK_review.RData",
     Data , nConds, nRatings,  maxrt, cond_levels,
     Data_RatingDist_corr_cond ,
     Data_MRating_corr_cond ,
     fits_WEVmodels ,  fits_RMmodels , fit_dynVib,fit_DDMConf,
     preds_RM ,  preds_WEV , pred_dynVib, pred_DDMConf, Ns_part,
     Preds_RatingDist_corr_cond ,  Preds_MRating_corr_cond ,
     Preds_RTdens_corr_cond_rating  ,
     Preds_RTQuants_corr_rating, Data_RTQuants_corr_rating,
     Preds_RTQuants_corr_cond , Data_RTQuants_corr_cond )

print("End of fitting and prediction. Saved all results.")

###############################################################



###### Comutation of marginalized log-likelihoods for specific subdomains
load("collected_fitsNpredicts_KonfMaskTime_review.RData")
predicted_resp_dist <- rbind(preds_WEV, preds_RM,
                             pred_dynVib[,names(preds_RM)],
                             pred_DDMConf[,names(preds_RM)])%>%
  mutate(rating=ifelse(model == "2DSD" & response==-1& participant ==10, ifelse(rating-1>0, rating-1, 5), rating)) %>%
  mutate(stimulus=ifelse(model %in% c("IRM", "IRMt", "PCRMt", "PCRM"), 2*stimulus-3, stimulus),
         response = ifelse(response=="lower", -1, ifelse(response=="upper", 1, response)),
         response = as.numeric(response),
         response=ifelse(model %in% c("IRM", "IRMt", "PCRMt", "PCRM"), 2*response-3, response))
stim_levels <- unique(Data$stimulus)
resp_dist_Data <- Data %>% group_by(participant, condition, stimulus, response, rating) %>%
  summarise(Ntrials = n()) %>%
  mutate(condition= as.numeric(condition), 
         stimulus = if_else(stimulus==stim_levels[1], 1, -1),
         response = if_else(response==stim_levels[1], 1, -1))
MargBIC <- left_join(resp_dist_Data, predicted_resp_dist) %>% 
  mutate(logprob = log(p)*Ntrials) %>%
  group_by(participant, model) %>%
  summarise(negloglik = -sum(logprob),
            Ntrials = sum(Ntrials),
            npar = case_when(model[1]=="DDMConf" ~ 19,
                             model[1]=="2DSD" ~ 20,
                             model[1]=="WEVmu" ~ 23,
                             model[1] %in% c("IRM", "PCRM") ~17,
                             TRUE ~ 19),
            BIC = 2*negloglik+log(Ntrials)*npar) %>%
  group_by(model) %>% 
  summarise(Mnegloglik = mean(negloglik), 
            SDnegloglik=sd(negloglik),
            MBIC=mean(BIC),
            SDBIC=sd(BIC))
MargBIC
