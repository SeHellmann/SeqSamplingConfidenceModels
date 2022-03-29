###########################################################################
#####      Model mimikry analysis for SeqSamplingConfModels         #######
###########################################################################

# Sebastian Hellmann, 21.10.2021

# 1) Load fits from file
# 2) For the best three models simulate data and fit those best models
# 3) Compare and visualize model mimikry


###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/AccumulatorModels")
# or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

{
  library(plyr)
  library(snow)
  library(doSNOW)
  library(BayesFactor)
  library(tidyverse)
  library(dynWEV)
  library(RColorBrewer)
  library(gridExtra)
}
#=============================================================
#=============================================================
###   Load results/fits to empirical data        ####
if (!file.exists("JRRDK/collected_fitsNpredicts.RData")) {
  stop("No file 'collected_fitsNpredicts.RData' in current folder;
       Maybe change working directory to right folder or run 'Script_SeqSampConfModels_....R' first!")
}
load("JRRDK/collected_fitsNpredicts.RData")
Data_JRRDK <- ungroup(Data)
Data_JRRDK$participant <- Data_JRRDK$participant + 100

#fits_WEVmodels <- fits_2DSD_WEV_neu
fits_WEVmodels_JRRDK <- fits_WEVmodels  %>% mutate(participant = participant +100)
fits_RMmodels_JRRDK <- fits_RMmodels %>% mutate(participant = participant +100)

load("KonfRDKtime/collected_fitsNpredicts.RData")
Data <- rbind(select(ungroup(Data), participant, condition, correct, rt, Rating, stimulus, response,  rating, nrows),
              select(Data_JRRDK, participant, condition, correct, rt, Rating, stimulus, response, rating, nrows))
fits_WEVmodels <- rbind(fits_WEVmodels, fits_WEVmodels_JRRDK)
fits_RMmodels <- rbind(fits_RMmodels, fits_RMmodels_JRRDK)

N_trials_part <- Data %>% group_by(participant) %>%
  summarise(Ntrials = n())

#=============================================================
#=============================================================
########### 1) Simulate and fit      ##############
# Define a function to be parallelized over for simulation and fitting
generative_models_for_comparison <- c("2DSD", "PCRMt")  # Choose models for comparison here
Nsims = 1
planned_analysis <- expand.grid(participant = unique(fits_RMmodels$participant),
                                n_simu = 1:Nsims,
                                model=generative_models_for_comparison) 

fun_sim_fit <- function(simu_row) {
  set.seed(simu_row$participant+100*simu_row$n_simu+20*str_length(simu_row$model))
  gen_model <- simu_row[['model']]
  cur_participant <- simu_row[['participant']]
  
  outnames <- union(colnames(fits_RMmodels), colnames(fits_WEVmodels))[-c(1,2)]
  out <- data.frame(matrix(NA, nrow=0, ncol=length(outnames)))
  colnames(out) <- outnames
  
  if (grepl(pattern = "RM", gen_model)) {
    paramDf <- subset(fits_RMmodels, participant == cur_participant & model==gen_model)
    sim_df <- rRM(paramDf = paramDf, model = gen_model,
                          n = ceiling(N_trials_part[N_trials_part$participant== cur_participant,][['Ntrials']]/(2*5)),     # number of trials per cond and stimulus  
                          delta=0.01,        # discretization step for simulation (in sec.) 
                          maxrt=30)         # maximum decision time for each trial
    sim_df <- sim_df %>% filter(response %in% c(1,2)) # if maxrt is exceeded; response is set to 0
    sim_df$participant <- cur_participant
    res_fit_gen_model <- fitRM(sim_df, model = gen_model, logging = TRUE, nRatings=5)
    neglogLik <- -LogLikRM(sim_df, paramDf = paramDf, model = gen_model)
  } else {
    paramDf <- subset(fits_WEVmodels, participant == cur_participant & model==gen_model)
    sim_df <- simulateWEV(paramDf = paramDf, model = gen_model, 
                          n = ceiling(N_trials_part[N_trials_part$participant== cur_participant,][['Ntrials']]/(2*5)),     # number of trials per cond and stimulus  
                          delta=0.01,        # discretization step for simulation (in sec.) 
                          maxrt=30)         # maximum decision time for each trial
    sim_df <- sim_df %>% mutate(rt = rt+ paramDf$tau)
    sim_df <- sim_df %>% filter(response %in% c(-1, 1)) # if maxrt is exceeded; response is set to 0
    simd_df <- subset(sim_df, response !=0)
    sim_df$participant <- cur_participant
    res_fit_gen_model <- fitWEV(sim_df, model = "2DSD", logging = TRUE, restr_tau="simult_conf", nRatings=5)
    neglogLik <- -LogLikWEV(sim_df, paramDf = paramDf, model = gen_model, simult_conf = TRUE)
  }
  out[1, colnames(res_fit_gen_model)] <- res_fit_gen_model
  
  res_fit_WEV <- fitWEV(sim_df, model = "WEVmu", logging = TRUE, restr_tau="simult_conf", nRatings=5)
  colnames(res_fit_WEV) <- paste("WEV", colnames(res_fit_WEV), sep="_")
  
  out <- cbind(out, res_fit_WEV)
  out$true_negLik <- neglogLik
  print(paste("Finished participant:", simu_row$participant, "model:", simu_row$model, "simulation:", simu_row$n_simu))
  return(out)
}


n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores, "SOCK", outfile = "")
registerDoSNOW(cl)
clusterEvalQ(cl, library(dynWEV))
clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("planned_analysis", "fun_sim_fit", "fits_WEVmodels", "fits_RMmodels", "N_trials_part"))


t00 <- Sys.time()
mimikry_collected_results <- ddply(planned_analysis,.(participant, model, n_simu),
                              .fun = fun_sim_fit,
                              .parallel = T)
# [1] "Mimikry analysis took... 29337.81  mins"
dir.create("autosave_mimikry")
save(file = "autosave_mimikry/results_mimikry_analysis.RData", mimikry_collected_results,fun_sim_fit)
print("Finished model mimikry analysis!")
print(paste("Mimikry analysis took...",
            as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
            " mins"))
stopCluster(cl)






########    Analyse Mimikry Results
load("autosave_mimikry/results_mimikry_analysis.RData")
head(mimikry_collected_results)


### How was the weight parameter fitted?
temp <- fits_WEVmodels %>% filter(model=="WEVmu") %>% 
  select(participant, w, BICWEV=BIC) %>% 
  left_join(select(filter(fits_WEVmodels, model=="2DSD"), participant, BIC)) 
ggplot(temp, aes(x=BIC, y=BICWEV/BIC, color=w))+
  geom_point(size=2.5) +
  geom_abline(slope=0, intercept = 1)+
  scale_color_continuous(low="red", high="blue")
min(fits_WEVmodels$w, na.rm = TRUE)
fits_WEVmodels[replace_na(fits_WEVmodels$w<0.4, FALSE),]
ggplot(temp, aes(x=w, y=BICWEV/BIC, color=w))+
  geom_point(size=2.5) +
  geom_abline(slope=0, intercept = 1)+
  scale_color_continuous(low="red", high="blue")
mimikry_collected_results %>% filter(model=="2DSD") %>%
  select(WEV_w) %>%
  summary(.)
mimikry_collected_results %>% filter(model=="2DSD") %>%
  select(WEV_w) %>%
  boxplot(.)

mimikry_collected_results %>% filter(model=="2DSD") %>%
  arrange(WEV_w)



mimikry_collected_results <- mimikry_collected_results %>%
  mutate(BIC_diff = WEV_BIC-BIC,
         true_likdiff = true_negLik - negLogLik) %>%
  mutate(Ntrials = as.character(factor(cut(N, breaks=c(0, 700,  1000, 1500, 3000)),
                                       labels=c(640, 880, 1280, 1920))))
mimikry_collected_results[, c("N", "Ntrials")]

labels_model <- c("Generative model: 2DSD", "Generative model: PCRMt")
names(labels_model) <- c("2DSD", "PCRMt")
BIC_diff_labels <- data.frame(x=0.6, y=c(3*log(640), 3*log(1920), 4*log(640), 4*log(1920)),
                              label=c("3*log(640)=19.38", "3*log(1920)=22.68", "4*log(640)=25.85", "4*log(1920)=30.24"),
                              color=c("640", "1920","640", "1920"),
                              model=c("2DSD", "2DSD", "PCRMt", "PCRMt"))
ggplot(mimikry_collected_results, aes(y=BIC_diff, x=model))+
  geom_violin()+
  geom_boxplot()+
  geom_jitter(aes(color=Ntrials),width=0.2, size=2)+
  geom_label(data=BIC_diff_labels, aes(x=x, y=y, label=label, color=color))+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_discrete(name="")+
  ggtitle("Results from model mimikry analysis for RDK experiments")+
  ylab(expression(BIC[dynWEV]-BIC[genmodel]))+
  facet_wrap(~model, scales="free", labeller = labeller(model=labels_model))+
  theme_bw()
