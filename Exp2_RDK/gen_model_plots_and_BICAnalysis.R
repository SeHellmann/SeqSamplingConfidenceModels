####################################################################
###  Generate figures for comparison between model predictions  #### 
###   and empirical data and conduct Bayesian tests for BIC     ####  
####################################################################
#### Structure:                                                 ####
#### 1) Load and combine data and model fits from two data sets ####
#### 2) Aggregate combined data and predictions for plots       ####
#### 3) Produce Figures for model fits                          ####
#### 4) BIC Analyses: Bayes t test + Figure for BIC differences ####
## Preamble:                                                    ####
rm(list = ls())
library(Rmisc)
library(ggpubr)
library(dynWEV)
library(tidyverse)
library(gridExtra)
library(grid)
library(BayesFactor)

## Include font style and define colors for plots   
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#1f78b4", "#a6cee3", "#33a02c","#b2df8a", "#e31a1c","#fb9a99")
two_colors_model <- c("#1f78b4", "#b2df8a")


####################################################################
#### 1) Load and combine data and model fits from two data sets ####
load("JRRDK/collected_fitsNpredicts.RData")
Data_JRRDK <- ungroup(Data)
Data_JRRDK$participant <- Data_JRRDK$participant + 100
preds_RM_JRRDK <- preds_RM %>% mutate(participant = participant +100)
preds_WEV_JRRDK <- preds_WEV %>% mutate(participant = participant +100)
RT_dist_JRRDK <- RT_dist %>% mutate(participant = participant +100)
fits_WEVmodels_JRRDK <- fits_WEVmodels  %>% mutate(participant = participant +100)
fits_RMmodels_JRRDK <- fits_RMmodels %>% mutate(participant = participant +100)

load("KonfRDKtime/collected_fitsNpredicts.RData")
Data <- rbind(select(ungroup(Data), participant, condition, correct, rt, Rating, stimulus, response,  rating, nrows),
              select(Data_JRRDK, participant, condition, correct, rt, Rating, stimulus, response, rating, nrows))
preds_RM <- rbind(preds_RM, preds_RM_JRRDK) 
preds_WEV <- rbind(preds_WEV, preds_WEV_JRRDK)
fits_WEVmodels <- rbind(fits_WEVmodels, fits_WEVmodels_JRRDK)
fits_RMmodels <- rbind(fits_RMmodels, fits_RMmodels_JRRDK)


### For discussion: Proportion of post-decisional accumulation in response time
Data %>% left_join(select(filter(fits_WEVmodels, model=="WEVmu"), participant, tau)) %>%
  summarise(prop_post = mean(tau/rt))




#### 2) Aggregate combined data and predictions for plots       ####
### (everything like in the script)       
#### Data                                                       ####
Data_RatingDist_part <- Data %>% 
  group_by(stimulus, condition, rating, response, correct, participant) %>% 
  summarise(p = n()/(mean(nrows))) %>% 
  full_join(y = expand.grid(stimulus = unique(Data$stimulus), 
                            condition = unique(Data$condition),
                            rating = 1:nRatings, 
                            correct = c(0,1), 
                            participant = unique(Data$participant))) %>%
  mutate(p = ifelse(is.na(p), 0, p)) %>%
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

#### Compute Mean Rating of the Data    
Data_MRating_corr_cond <- Data_RatingDist_part %>% 
  group_by(condition, participant, correct) %>%
  summarise(MRating = sum(rating*p)/sum(p)) %>%
  group_by(condition, correct) %>%
  summarise(SER = sd(MRating, na.rm = TRUE)/sqrt(n()),
            MRating = mean(MRating, na.rm = TRUE))

#### Compute Reaction Time Quantiles of the Data for different accuracy and ratings 
# We do not respect other factors, here, but handle the data set as one distribution and
# only split between incorrect and correct decisions and confidence ratings.
Data_RTQuants_corr_rating <- Data %>%
  group_by(rating, correct) %>%
  summarise(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 


#### Predictions                                                ####
##   Combine and Aggregate Confidence Rating Distribution   ##
Preds_RatingDist_corr_cond <- rbind(preds_WEV, preds_RM) %>%
  left_join(Ns_part) %>%
  group_by(model, rating, correct, condition) %>%
  summarise(p = mean(p)) %>%
  mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
                        labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")), 
         condition = factor(condition,levels = 1:nConds, 
                            labels =  as.character(cond_levels*100)))

### Sanity checks:
# Preds_RatingDist_corr_cond %>% filter(model %in% c("IRM", "IRMt", "PCRM", "PCRMt")) %>%
#   group_by(model, condition) %>% summarise(p = sum(p))
# Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
#   summarise(p=sum(p))

###    Compute Mean Rating Accross Conditions 
# This is good to visualize a folded-X- and double-increase-pattern 
Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% group_by(model, condition, correct) %>%
  summarise(MRating = sum(p*rating)/sum(p))


####### Compute predicted RT-distribution on a common interval  ####
## As predictions of the density are made in discrete steps, 
## it is easier to have the same grid for all models and 
## participants. So, we compute the density prediction again 
## for the combined data to have a common grid. 
## 
## EITHER: Load previously computed density (saves a lot of time!) ####
## As the computations take some time, we load the results 
## form the computations. The code is commented below and 
## may be run if necessary. 
maxrt <- max(Data$rt)  
try(load("common_RDK_RTpreds.RData"))

## OR: Recomputation of predicted RT density for common grid (takes some time!) ####
# compute_RTdens_from_fits <- function(Confpred_data) {
#   if (Confpred_data$model[1] %in% c("IRM","PCRM","IRMt", "PCRMt")) {
#     paramDf <- subset(fits_RMmodels, 
#                       model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
#     res <- predictRM_RT(paramDf = paramDf, 
#                         model = Confpred_data$model[1], 
#                         maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0),
#                         scaled = TRUE, DistConf = Confpred_data,
#                         .progress = FALSE)
#   } else {
#     paramDf <- subset(fits_WEVmodels, 
#                       model == Confpred_data$model[1] & participant == Confpred_data$participant[1])
#     res <- predictWEV_RT(paramDf = paramDf, 
#                          model = Confpred_data$model[1], 
#                          maxrt = maxrt, subdivisions = 300, minrt = min(fits_WEVmodels$t0, fits_RMmodels$t0), 
#                          scaled = TRUE, DistConf = Confpred_data, simult_conf = TRUE,
#                          .progress = FALSE)
#   }
#   return(res)
# }
# n.cores <- parallel::detectCores() - 1
# cl <- makeCluster(n.cores, "SOCK", outfile = "")
# registerDoSNOW(cl)
# clusterEvalQ(cl, library(dynWEV))
# clusterExport(cl, c("fits_WEVmodels"))
# clusterExport(cl, c("fits_RMmodels"))
# clusterExport(cl, c("preds_WEV", "preds_RM", "maxrt", "compute_RTdens_from_fits"))
# RT_dist_RM <- ddply(preds_RM,.(participant, model),
#                     .fun = compute_RTdens_from_fits,
#                     .parallel = T)
# RT_dist_RM <- RT_dist_RM %>% group_by(model, participant, correct, rating, condition, rt) %>%
#   summarise(dens = mean(dens), densscaled = mean(densscaled))
# 
# RT_dist_WEV <- ddply(preds_WEV,.(participant, model),
#                      .fun = compute_RTdens_from_fits, 
#                      .parallel = T)
# RT_dist_WEV <- RT_dist_WEV %>% group_by(model, participant, correct, rating, condition, rt) %>%
#   summarise(dens = mean(dens), densscaled = mean(densscaled))
# RT_dist <- rbind(RT_dist_WEV, RT_dist_RM)
# stopCluster(cl)
# rm(list = c("cl"))
# 
# ####   Combine and Aggregate RT densities over participants
# ## For an aggregation, which is equivalent to the aggregation of 
# ## empirical observations, we compute a weighted mean of participant 
# ## densities with trial numbers as weights. 
# Ns_part <- Data %>% group_by(participant) %>% summarise(N=n(), MinRT = min(rt))  %>%
#   select(participant, N)
# Preds_RTdens_corr_cond_rating <-  RT_dist %>%
#   left_join(Ns_part) %>%
#   ungroup() %>%
#   select(-participant) %>%
#   group_by(rating, condition, model, correct, rt) %>%
#   summarise(dens = mean(dens), densscaled = mean(densscaled)) %>%
#   mutate(model = factor(model, levels = c("WEVmu", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM"),
#                         labels = c("dynWEV", "2DSD", "IRMt", "IRM", "PCRMt", "PCRM")),
#          condition = factor(condition,levels = 1:nConds,
#                             labels = as.character(cond_levels*100)))
# 
# save(file="common_RDK_RTpreds.RData",
#      Preds_RTdens_corr_cond_rating, Ns_part)


####### Compute aggregated distribution summaries for RT-dist   ####

####    Computation of RT quantiles seperated by accuracy and confidence rating
Preds_RTQuants_corr_rating <- Preds_RTdens_corr_cond_rating %>% 
  group_by(model, rt, correct, rating) %>%
  summarise(dens = mean(dens))%>%
  group_by(model, correct, rating) %>% 
  do(RTDensityToQuantiles(.)) 

####### Save collected data and predictions for both experiments####
save(file="collected_fitsNpredicts_bothRDKExperiments.RData", 
     Data , nConds, nRatings,  maxrt, cond_levels,
     Data_RatingDist_corr_cond,  Data_MRating_corr_cond,   
     fits_WEVmodels ,  fits_RMmodels ,
     preds_RM ,  preds_WEV ,
     Preds_RatingDist_corr_cond ,  Preds_MRating_corr_cond ,  
     Preds_RTdens_corr_cond_rating, 
     Preds_RTQuants_corr_rating, 
     Data_RTQuants_corr_rating)



#### 3) Produce Figures for model fits                          ####
####### Figure 8: Plot of Mean Ratings                          ####
levels(Preds_MRating_corr_cond$model) <- c("dynamical weighted evidence\nand visibility", "two-stage dynamical\nsignal detection",
                                           "independent race model\n(time-dependent confidence)","independent race model\n(Balance of Evidence)",
                                           "anti-correlated race model\n(time-dependent confidence)","anti-correlated race model\n(Balance of Evidence)")
Data_MRating_corr_cond <- Data %>% group_by(participant, condition, correct) %>%
  summarise(MRating = mean(rating)) %>%
  summarySEwithin(measurevar = "MRating", idvar = "participant", withinvars = c("condition", "correct")) %>%
  rename(sewithin = se) %>% select(-MRating) %>%
  right_join(mutate(Data_MRating_corr_cond, correct=as.factor(correct)), by=c("condition", "correct"))

pd <- position_dodge(0.02) 
p_MRating <- ggplot(mutate(Data_MRating_corr_cond, condition = as.integer(condition)),
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  geom_line(data=Preds_MRating_corr_cond, aes(color=as.factor(correct)), size=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.2,position =pd, size=0.4) +
  geom_point(position = pd, fill="white")+
  facet_wrap(.~model, nrow = 3)+ #dir="v"
  ylab("Mean Confidence")+
  scale_x_discrete(name="Motion Coherence [%]")+   
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted",
                     labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed",
                     labels=c("Correct", "Wrong"))  +
  scale_y_continuous(limits = c(1, 5), breaks=1:5, labels = c("0", "25", "50", "75", "100"), name="Mean Confidence")+
  theme_bw() + 
  guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom",legend.direction = "vertical", #c(.005, .995), legend.justification = c(0,1),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
p_MRating
ggsave("figures/meanRating2.eps", 
       width = 17.62, height=14, units="cm",dpi=1200, device = cairo_ps)


####### Suppl Figure 2: Plot of Ratings Dist.                   ####
ann_text <- expand.grid(correct=c(0,1),
                        condition= levels(Preds_RatingDist_corr_cond$condition ),
                        label=NA, rating=0, p=1)
ann_text$label <- c("A", rep("", 9))
p_ratingdist <- ggplot(data=Data_RatingDist_corr_cond, aes(x=rating, y=p))+
  geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
  geom_point(data=Preds_RatingDist_corr_cond,aes(shape=model, fill=model), position=position_dodge(1), col="black")+
  facet_grid(cols=vars(correct), rows = vars(condition), 
             labeller = labeller(correct=c("0"="Wrong", "1"="Correct"), 
                                 condition=function(x) paste(x, "%")))+
  scale_x_continuous(
    name = "Identification confidence [% scale width]", 
    breaks = 1:5, 
    labels = c("0-20","20-40", "40-60", "60-80", "80-100")) +
  scale_y_continuous(name = "Probability", breaks = c(0, 0.4, 0.8))+ 
  #ggtitle("Confidence Distribution: Data vs. Predictions; Fitted per participants (aggregated for plotting)") +
  # scale_fill_brewer(type = "qual", palette=6, 
  #                   name = "Model prediction" ) +
  scale_shape_manual(name = "Model prediction",values=c(21, 22, 23, 24, 25, 21)) +
  scale_fill_manual(name = "Model prediction",values=model_colors) +
  theme_bw() + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0, size=9, family = "Times"),
        panel.spacing=unit(0, "lines"))+
  coord_cartesian(xlim = c(0.5,5.5),ylim=c(0, 0.9),  clip = 'off')
  # geom_text(data=ann_text, label=ann_text$label,size=14/.pt, 
  #          hjust=1, vjust=0, family="Times", fontface="bold")
p_ratingdist
ggsave("figures/RespDist2.eps", 
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)



####### Figure 9: RTQuantiles accross correct X rating          ####
Data_RTQuants_corr_rating <- Data_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_rating <- Preds_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
levels(Preds_RTQuants_corr_rating$model) <- c("dynamical weighted evidence\nand visibility", "two-stage dynamical\nsignal detection",
                                              "independent race model\n(time-dependent confidence)","independent race model\n(Balance of Evidence)",
                                              "anti-correlated race model\n(time-dependent confidence)","anti-correlated race model\n(Balance of Evidence)")

ggplot()+
  geom_point( data=Preds_RTQuants_corr_rating, aes(x=X, y=q,shape="Model fit"),
              size=1.5, fill="white", color="white") +  
  geom_line(data=Data_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Observed", color="Observed"), size=0.7)+
  geom_line(data=Preds_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Model fit", color="Model fit"), size=0.7)+

  geom_point(data=Data_RTQuants_corr_rating, aes(x=X, y=q, shape="Observed"),
             size=1.2)+
  scale_shape_manual(name="class", breaks = c("Model fit", "Observed"), values=c(32, 17))+ 
  scale_x_discrete(name="Confidence x Accuracy", 
                   labels=c("100-80", "60-80","40-60\n Wrong", "20-40",  "0-20", 
                            "0-20","20-40", "40-60\n Correct","60-80", "80-100"))+
  scale_y_continuous(breaks = c(2,4,6,8, 10), name="Reaction Time Quantiles [s]")+
  scale_linetype_manual(name="class", values= c("solid", "dashed"), breaks = c("Model fit", "Observed"))+
  scale_color_manual(name="class", values= two_colors_model, breaks = c("Model fit", "Observed"))+
  facet_wrap(.~model, nrow=3)+ #,dir="v"
  theme_bw() + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        legend.text = element_text(size=9),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(size=7, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9, family = "Times"),
        #legend.box = "horizontal",
        legend.position = "bottom", #c(.998, .998), legend.justification = c(1,1),#legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"),
        panel.spacing=unit(0, "lines"))
ggsave("figures/RTQuantsConf2.eps", 
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ungroup(Data) %>% filter(rating==5 & correct==0) %>% summarise(N=n()/nrow(Data))






#### 4) BIC Analyses: Bayes t test + Figure for BIC differences ####

# See which model performed best for how many participants
BIC_ranks <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                   fits_WEVmodels[,c("participant", "model", "BIC")])  %>% 
  arrange(participant, BIC) %>%
  mutate(rank = rep(1:6, length(unique(participant))))
print("Model ranks in terms of BIC: ")
BIC_ranks %>% group_by(model, rank) %>%
  summarise(N=n())

## Compute BIC differences between dynWEV and other models
descr_BIC <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                   fits_WEVmodels[,c("participant", "model", "BIC")]) %>%
  filter(model!="WEVmu") %>% 
  left_join(filter(fits_WEVmodels[,c("participant","model", "BIC")], model=="WEVmu"), by="participant") %>%
  group_by(model.x) %>% 
  summarise(MeanBIC = mean(BIC.x-BIC.y),
            SDBIC = sd(BIC.x-BIC.y))
descr_BIC

####### Bayes t-test between dynWEV and all other models        ####
compute_bf <- function(x,y) {
  drop <- is.na(x) | is.na(y)
  x <- x[!drop]
  y <- y[!drop]
  bf <- ttestBF(x,y,  paired = TRUE, rscale=1)
  res <-  data.frame(bf = as.numeric(extractBF(bf)["bf"]), Lower = NA, Upper = NA)
  if (! is.na(res$bf)) {
    posterior <- quantile(posterior(bf, iterations=100000)[,"delta"], probs = c(.025, .975))
    res$Lower <- posterior[1]
    res$Upper <- posterior[2]
  }
  res
}

criteria_comparison <- rbind(select(fits_WEVmodels, c("participant","model","BIC", "AICc", "AIC")),
                             select(fits_RMmodels, c("participant","model","BIC", "AICc", "AIC"))) %>% 
  pivot_longer(cols=c("BIC", "AICc", "AIC"), names_to="criteria", values_to="value") %>%
  #pivot_wider(id_cols = c("participant", "criteria"), names_from = model, values_from = c("value")) %>%
  filter(model != "WEVmu") %>%
  group_by(criteria, model) %>%
  summarise(compute_bf(x=value, 
                       y=subset(fits_WEVmodels, model =="WEVmu")[,cur_data_all()$criteria[1]]),
            M_diff=mean(value-subset(fits_WEVmodels, model =="WEVmu")[,cur_data_all()$criteria[1]], na.rm=TRUE)) %>%
  rename(BF_model_vs_WEVmu = bf)
criteria_comparison <- arrange(criteria_comparison, model, criteria)
criteria_comparison
knitr::kable(criteria_comparison, digits=3)



####### SupplFig 4: Binned BIC diffs btw dynWEV and other models#### 
pal = rev(c(rgb(50,3,10,maxColorValue = 255), rgb(128,34,45,maxColorValue = 255),
            rgb(180,68,78,maxColorValue = 255), 
            rgb(252,252,252,maxColorValue = 255),
            rgb(189,209,230,maxColorValue = 255),rgb(143,169,200,maxColorValue = 255),
            rgb(67, 147, 195,maxColorValue = 255)))
pal <- c(hsv(245/360, 100/100, 54/100), hsv(245/360, 67/100, 60/100),
         hsv(245/360, 51/100, 96/100),
         hsv(0, 0, 1),
         hsv(0,38/100, 100/100), hsv(0,67/100, 92/100), 
         hsv(0,100/100, 72/100))

plotBIC_data <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                      fits_WEVmodels[,c("participant", "model", "BIC")]) %>%
  filter(model !="WEVmu") %>%
  left_join(rename(fits_WEVmodels[fits_WEVmodels$model=="WEVmu" ,c("participant", "BIC")], BICWEV = BIC)) %>%
  mutate(model= factor(model, levels=c("PCRMt", "PCRM", "IRMt", "IRM",  "2DSD"),
                       labels=c("anti-correlated race model\n(time-dependent confidence)",
                                "anti-correlated race model\n(Balance of Evidence)",
                                "independent race model\n(time-dependent confidence)",
                                "independent race model\n(Balance of Evidence)",
                                "two-stage dynamical\nsignal detection")),
         deltaBIC = BIC-BICWEV,
         binnedDeltaBIC = as.integer(cut(deltaBIC, breaks=c(-Inf,-100, -10,-2,2, 10, 100, Inf))))
table(plotBIC_data$binnedDeltaBIC, plotBIC_data$model)
plotBIC_data <- plotBIC_data %>% left_join(summarise(group_by(Data, participant), N = n()))
# ggplot(plotBIC_data) +
#   geom_jitter(aes(x=N, y=deltaBIC, col=model))



Plot_deltaBIC_Confidence_2 <- 
  ggplot(plotBIC_data, 
         aes(x=factor(binnedDeltaBIC, levels=1:6, labels = 1:6), 
             fill= factor(binnedDeltaBIC, levels=1:7, labels = 1:7))) + 
  facet_wrap( ~ model, ncol=3, dir = "v",  scales = "fixed", as.table = F) + 
  geom_bar(colour = "black") + 
  scale_y_continuous(name = "Number of participants") + 
  scale_x_discrete(drop=FALSE, 
                   name = expression(Delta~BIC~vs.~dynamical~evidence~and~visibility~model))  + 
  scale_fill_manual(
    values=pal, 
    labels = c(expression(paste(Delta~BIC, " < ", "\u2212",100)),
               expression(paste("\u2212",100," < ", Delta~BIC, " < ", "\u2212",10)), 
               expression(paste("\u2212",10," < ", Delta~BIC, " < ", "\u2212",2)), 
               expression(paste("\u2212",2," < ", Delta~BIC, " < ", 2)),
               expression(paste(2," < ", Delta~BIC, " < ", 10)), 
               expression(paste(10, " < ", Delta~BIC, " < ", 100)), 
               expression(paste(Delta~BIC, " > ", 100))), 
    name = expression(Delta~BIC), drop=FALSE) + 
  theme_bw() + 
  theme(text = element_text(size=9, family="Times" ),
        strip.text = element_text(size=9, family = "Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text( angle=90),# switch off minor gridlines
        legend.text.align = 0, legend.text = element_text(size=9),
        panel.grid.minor = element_blank(), 
        plot.margin = margin(0, 0, 0, 0, "cm"),#  plot.margin = unit(c(.5,.5,.5,.5), "lines"),
        panel.spacing = unit(.0, "lines"), 
        legend.position = c(.125, .25), 
        legend.title = element_blank(),
        legend.key.width=unit(1,"line"),
        legend.key.height = unit(1, "line"))
# legend.title = element_blank(), 
# legend.key = element_blank(),
# legend.key.width=unit(2.5,"line"),
Plot_deltaBIC_Confidence_2


# dir.create("figures")
ggsave(file = "figures/BICs2_fewcats.png",
       dpi=600,
       height = 10, width = 17.62,
       units = "cm")
ggsave("figures/BICs2.eps", 
       width=17.62, height=10, dpi=1200, units="cm", device=cairo_ps)