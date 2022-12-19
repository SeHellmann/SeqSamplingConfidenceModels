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
library(Hmisc)
library(ggpubr)
library(dynWEV)
library(tidyverse)
library(gridExtra)
library(grid)
library(BayesFactor)
library(xtable)

## Include font style and define colors for plots
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#1f78b4", "#a6cee3", "#33a02c","#b2df8a", "#e31a1c","#fb9a99", "#000000", "#000000")
two_colors_model <- c("#1f78b4", "#b2df8a")


#source("CombineDataandFits_bothRDKs.R")

load("../Additional_Analyses/collected_fitsNpredicts_RDK_review.RData")
#### 3) Produce Figures for model fits                          ####

####### Supplement Figure: Plot of Fitted Accuracy              ####
Data_Acc <- Data %>% group_by(participant, condition) %>%
  summarise(Acc = mean(correct)) %>%
  summarySEwithin(measurevar = "Acc", idvar = "participant", withinvars = c("condition")) %>%
  rename(sewithin = se)
Preds_Acc <- Preds_RatingDist_corr_cond %>% group_by(model, condition) %>%
  summarise(Acc = sum(p*correct)/sum(p)) %>%
  mutate(condition = factor(condition, labels= levels(Data_Acc$condition)),
         model= factor(model, levels=c("dynWEV", "2DSD", "dynVib", "IRMt", "IRM", "DDMConf", "PCRMt", "PCRM"),
                       labels=c("dynamical weighted evidence\nand visibility (dynWEV)",
                                "two-stage dynamical\nsignal detection (2DSD)",
                                "dynamical visibility\n(dynVis)",
                                "independent race model\n(time-dependent confidence, IRMt)","independent race model\n(Balance of Evidence, IRM)",
                                "drift diffusion confidence\n(DDConf)",
                                "anti-correlated race model\n(time-dependent confidence, PCRMt)","anti-correlated race model\n(Balance of Evidence, PCRM)")))
p_Acc <- ggplot(Data_Acc,
                aes(x=condition, y=Acc)) +
  geom_line(data=Preds_Acc, aes(x = as.numeric(condition), linetype="Predicted"), size=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=Acc-sewithin, ymax=Acc+sewithin), colour="black", width=.2, size=0.4) +
  geom_point(fill="white", aes(shape="Observed"))+
  facet_wrap(.~model, nrow = 3)+ #, dir="v"
  ylab("Mean accuracy")+
  scale_x_discrete(name="Motion coherence [%]")+
  scale_linetype_manual(name="", values=1) +
  scale_shape_manual(values=c(21),name = "")  +
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = c(.80, .005), legend.justification = c(0,0),
        legend.margin = margin(0, 0, 0, 0, "cm"),#legend.direction = "horicontal",
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
p_Acc
dir.create("figures", showWarnings = FALSE)
ggsave("figures/modelAccuracy2.eps",
       width = 17.62, height=14, units="cm",dpi=1200, device = cairo_ps)


####### Figure 8: Plot of Mean Ratings                          ####
Preds_MRating_corr_cond$model <- factor(Preds_MRating_corr_cond$model, 
                                        levels=c("dynWEV", "2DSD", "dynVib", "IRMt", "IRM", "DDMConf", "PCRMt", "PCRM"),
                                        labels=c("dynamical weighted evidence\nand visibility (dynWEV)",
                                                 "two-stage dynamical\nsignal detection (2DSD)",
                                                 "dynamical visibility\n(dynVis)",
                                                 "independent race model\n(time-dependent confidence, IRMt)","independent race model\n(Balance of Evidence, IRM)",
                                                 "drift diffusion confidence\n(DDConf)",
                                                 "anti-correlated race model\n(time-dependent confidence, PCRMt)","anti-correlated race model\n(Balance of Evidence, PCRM)"))

Data_MRating_corr_cond <- Data %>% group_by(participant, condition, correct) %>%
  summarise(MRating = mean(rating)) %>%
  summarySEwithin(measurevar = "MRating", idvar = "participant", withinvars = c("condition", "correct")) %>%
  rename(sewithin = se) %>% select(-MRating) %>%
  right_join(mutate(Data_MRating_corr_cond, correct=as.factor(correct)), by=c("condition", "correct"))

pd <- position_dodge(0.02)
p_MRating <- ggplot(mutate(Data_MRating_corr_cond, condition = as.integer(condition)),
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  geom_ribbon(data=Preds_MRating_corr_cond, 
              aes(ymin=MRating-sewithin, ymax=MRating+sewithin, fill=as.factor(correct)), alpha=0.5)+
  geom_line(data=Preds_MRating_corr_cond, aes(color=as.factor(correct)), size=0.8)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.25,position =pd, size=0.6) +
  geom_point(position = pd, fill="white", size=1.8)+
  facet_wrap(.~model, nrow = 3)+ #dir="v"
  scale_x_discrete(name="Motion coherence [%]")+
  scale_fill_manual(values= two_colors_correct, breaks=c(1,0),
                    name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted",
                     labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed",
                     labels=c("Correct", "Wrong"))  +
  scale_y_continuous(limits = c(1, 5), breaks=1:5, labels = c("0", "25", "50", "75", "100"), name="Mean confidence")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3), fill="none")+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = c(.70, .005), legend.justification = c(0,0),
        legend.direction = "vertical", #legend.position = "bottom"
        legend.box="horizontal",
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.key.width=unit(1,"line"),
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
    name = "Confidence [% scale width]",
    breaks = 1:5,
    labels = c("0-20","20-40", "40-60", "60-80", "80-100")) +
  scale_y_continuous(name = "Probability", breaks = c(0, 0.4, 0.8))+
  #ggtitle("Confidence Distribution: Data vs. Predictions; Fitted per participants (aggregated for plotting)") +
  # scale_fill_brewer(type = "qual", palette=6,
  #                   name = "Model prediction" ) +
  scale_shape_manual(name = "Model prediction",values=c(21, 22, 23, 24, 25, 21, 21, 24)) +
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
       width=17.62, height=20, dpi=1200, units="cm", device=cairo_ps)



####### Figure 9: RTQuantiles accross correct X rating          ####
Data_RTQuants_corr_rating <- Data_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_rating <- Preds_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_rating<- ungroup(Preds_RTQuants_corr_rating)
Data_RTQuants_corr_rating <- ungroup(Data_RTQuants_corr_rating)
Preds_RTQuants_corr_rating$model1 <- factor(Preds_RTQuants_corr_rating$model, 
                                            levels=c("dynWEV", "2DSD", "dynVib", "DDMConf", "IRMt", "IRM", "PCRMt", "PCRM"),
                                            labels=c("dynamical weighted\nevidence and\nvisibility (dynWEV)",
                                                     "two-stage dynamical\nsignal detection\n(2DSD)",
                                                     "dynamical visibility\n(dynVis)",
                                                     "drift diffusion\nconfidence (DDConf)",
                                                     "independent race\nmodel (time-dependent\nconfidence, IRMt)",
                                                     "independent race\nmodel (Balance of\nEvidence, IRM)",
                                                     "anti-correlated race\nmodel (time-dependent\nconfidence, PCRMt)",
                                                     "anti-correlated race\nmodel (Balance of\nEvidence, PCRM)"))
ggplot()+
  geom_line(data=mutate(Preds_RTQuants_corr_rating, correct=factor(correct, labels=c("Wrong", "Correct")),
                        rating = as.factor(rating)),
            aes(x=rating, y=log(q), group=as.factor(p),color=correct), size=0.7)+
  geom_point(data=mutate(Data_RTQuants_corr_rating, correct=factor(correct, labels=c("Wrong", "Correct")),
                         rating = as.factor(rating)),
             aes(x=rating, y=log(q), shape=correct), size=1.2, fill="white")+
  scale_color_manual(values= two_colors_correct, breaks=c("Correct", "Wrong"),
                     name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c("Correct", "Wrong"),
                     name = "Observed", labels=c("Correct", "Wrong"))  +
  scale_x_discrete(name="Confidence [% scale width]", breaks=1:5,
                   labels=c("0-20","20-40", "40-60","60-80", "80-100"))+
  scale_y_continuous(breaks = log(c(2, 3, 5, 8)), 
                     labels = paste(c(2, 3, 5, 8), sep=""), 
                     name="Reaction Time Quantiles [s] (log scaled)")+
  facet_grid(cols=vars(correct), rows=vars(model1),
             labeller = )+ #,dir="v"
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        legend.text = element_text(size=9),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(size=7, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9, family = "Times", angle=0),
        strip.text.y.right = element_text(angle = 0),
        legend.box = "horizontal",
        legend.position = "bottom",# legend.position = c(.70, .005), legend.justification = c(0,0),
        legend.direction = "horizontal", legend.spacing.y = unit(0, "lines"),
        legend.margin =margin(0,0,0,1, "cm"), legend.box.spacing = unit(0.2,"lines"),
        #legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.width=unit(1.5,"line"),
        panel.spacing=unit(0, "lines"))
ggsave("figures/RTQuantsConf2.eps",
       width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
ggsave("../../Draft/figures/results/RTQuantsConf2.eps",
       width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
ggsave("../../Draft/figures/results/RTQuantsConf2.eps",
       width=18.288, height=15.24, dpi=1200, units="cm", device=cairo_ps)
ungroup(Data) %>% filter(rating==5 & correct==0) %>% summarise(N=n()/nrow(Data))






#### 4) BIC Analyses: Bayes t test + Figure for BIC differences ####

# See which model performed best for how many participants
BIC_ranks <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                   fits_WEVmodels[,c("participant", "model", "BIC")],
                   fit_dynVib[,c("participant", "model", "BIC")],
                   fit_DDMConf[,c("participant", "model", "BIC")])  %>%
  arrange(participant, BIC) %>%
  mutate(rank = rep(1:8, length(unique(participant))))
print("Model ranks in terms of BIC: ")
summary_BIC_ranks <- BIC_ranks %>% group_by(model, rank) %>%
  summarise(N=n()) 

## Compute BIC differences between dynWEV and other models
descr_BIC <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                   fits_WEVmodels[,c("participant", "model", "BIC")],
                   fit_dynVib[,c("participant", "model", "BIC")],
                   fit_DDMConf[,c("participant", "model", "BIC")]) %>%
  filter(model!="WEVmu") %>%
  left_join(filter(fits_WEVmodels[,c("participant","model", "BIC")], model=="WEVmu"), by="participant") %>%
  group_by(model.x) %>%
  summarise(MeanBIC = mean(BIC.x-BIC.y, na.rm=TRUE),
            SDBIC = sd(BIC.x-BIC.y, na.rm=TRUE))
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
                             select(fits_RMmodels, c("participant","model","BIC", "AICc", "AIC")),
                             arrange(select(fit_dynVib, c("participant","model","BIC", "AICc", "AIC")),participant),
                             arrange(select(fit_DDMConf, c("participant","model","BIC", "AICc", "AIC")),participant)) %>%
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
                      fits_WEVmodels[,c("participant", "model", "BIC")],
                      arrange(select(fit_dynVib, c("participant","model","BIC")),participant),
                      arrange(select(fit_DDMConf, c("participant","model","BIC")),participant)) %>%
  filter(model !="WEVmu") %>%
  left_join(rename(fits_WEVmodels[fits_WEVmodels$model=="WEVmu" ,c("participant", "BIC")], BICWEV = BIC)) %>%
  mutate(model= factor(model, levels=c("DDMConf", "IRMt", "PCRMt", "2DSD", "IRM",  "PCRM",  "dynVib"),
                       labels=c("drift diffusion\nconfidence (DDConf)",
                                "independent race model (time-\ndependent confidence, IRMt)",
                                "anti-correlated race model (time-\ndependent confidence, PCRMt)",
                                "two-stage dynamical\nsignal detection (2DSD)",
                                "independent race model\n(Balance of Evidence, IRM)",
                                "anti-correlated race model\n(Balance of Evidence, PCRM)",
                                "dynamical visibility\n(dynVis)")),
         deltaBIC = BIC-BICWEV,
         binnedDeltaBIC = as.integer(cut(deltaBIC, breaks=c(-Inf,-100, -10,-2,2, 10, 100, Inf))))
table(plotBIC_data$binnedDeltaBIC, plotBIC_data$model)
plotBIC_data <- plotBIC_data %>% left_join(summarise(group_by(Data, participant), N = n()))
# ggplot(plotBIC_data) +
#   geom_jitter(aes(x=N, y=deltaBIC, col=model))


Plot_deltaBIC_Confidence_2 <-
  ggplot(subset(plotBIC_data,!is.na(binnedDeltaBIC)),
         aes(x=factor(binnedDeltaBIC, levels=1:6, labels = 1:6),
             fill= factor(binnedDeltaBIC, levels=1:7, labels = 1:7))) +
  facet_wrap( ~ model, ncol=3, dir = "v",  scales = "fixed", as.table = T) +
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
        legend.position = c(.825, .25),
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


### BIC analysis for two experimental periods separately
fits_WEVmodels$exp_name <- if_else(fits_WEVmodels$participant <100, "KonfRDKtime", "JRRDK")
criteria_comparison_perperiod <- rbind(select(fits_WEVmodels, c("participant","model","BIC", "AICc", "AIC")),
                             select(fits_RMmodels, c("participant","model","BIC", "AICc", "AIC")),
                             arrange(select(fit_dynVib, c("participant","model","BIC", "AICc", "AIC")),participant),
                             arrange(select(fit_DDMConf, c("participant","model","BIC", "AICc", "AIC")),participant)) %>%
  pivot_longer(cols=c("BIC", "AICc", "AIC"), names_to="criteria", values_to="value") %>%
  #pivot_wider(id_cols = c("participant", "criteria"), names_from = model, values_from = c("value")) %>%
  filter(model != "WEVmu") %>%
  mutate(exp_name = if_else(participant <100, "KonfRDKtime", "JRRDK")) %>%
  group_by(criteria, model, exp_name) %>%
  summarise(compute_bf(x=value,
                       y=subset(fits_WEVmodels, model =="WEVmu" & exp_name==cur_data_all()$exp_name[1])[,cur_data_all()$criteria[1]]),
            M_diff=mean(value-subset(fits_WEVmodels, model =="WEVmu"& exp_name==cur_data_all()$exp_name[1])[,cur_data_all()$criteria[1]], na.rm=TRUE)) %>%
  rename(BF_model_vs_WEVmu = bf)
criteria_comparison_perperiod <- arrange(criteria_comparison_perperiod, model, criteria)
criteria_comparison_perperiod
knitr::kable(arrange(subset(criteria_comparison_perperiod, criteria == "BIC"),exp_name), digits=3)













####### Produce Tables for Supplementary Material      #####

### Table with mean Parameter estimations

Fits <- bind_rows(fits_WEVmodels, rename(fits_RMmodels, A=a, B=b),
                  fit_DDMConf, fit_dynVib)
meanFits <- Fits %>%
  select(-c("participant", "k", "N", "AICc", "negLogLik", "fixed")) %>%
  group_by(model) %>%
  summarise(across(.cols = everything(),
                   ~paste(format(round(mean(.x[abs(.x)<1e+24], na.rm=TRUE),2),nsmall=2),
                          " (", format(round(sd(.x[abs(.x) <1e+24], na.rm = TRUE), 2), nsmall=2), ")", sep=""))) %>%
  select(c("model", paste0("v", 1:5), "sv", "a", "z", "sz", "A", "B", paste0("thetaUpper", 1:4), paste0("thetaLower", 1:4),
           "t0", "st0", "tau", "w","sig", "sigmu", "wx", "wrt", "wint"))
names(meanFits)[2:6] <- paste0("$\\nu_", 1:5, "$")
names(meanFits)[7] <- paste0("$s\\nu$")
names(meanFits)[10] <- paste0("$s_z$")
names(meanFits)[13:16] <- paste0("$\\theta_{1,", 1:4, "}$")
names(meanFits)[17:20] <- paste0("$\\theta_{-1,", 1:4, "}$")
names(meanFits)[21:22] <- c("$t_0$", "$s_{t0}$")
names(meanFits)[23] <- paste0("$\\tau$")
names(meanFits)[25: 26] <- c("$s_V$", "$\\sigma_V$")
names(meanFits)[27:29] <- paste0("$w_{", c("x", "rt", "int"), "}$")
names(meanFits)
meanFits[meanFits=="NA (NaN)"]
meanFits[meanFits=="NA (NA)"]
meanFits[meanFits=="NaN (NaN)"]
meanFits[meanFits=="NaN (NA)"] <- "--"
meanFits[meanFits$model %in% c("IRM", "PCRM"),27:29] <- "--"
meanFits[meanFits$model %in% c("dynVib"),23:24] <- "--"



### Use this to produce a wide table
# meanFits <- column_to_rownames(meanFits, "model")
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="../Supplement/figures/tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)

### Use this to produce a long table
meanFits_long <- meanFits %>% pivot_longer(cols = 2:ncol(meanFits), names_to = "Parameter") %>%
  pivot_wider(id_cols = "Parameter", names_from = "model", values_from = value) %>%
  select(Parameter, dynWEV = WEVmu, '2DSD', DDConf=DDMConf, dynVis=dynVib, IRM, IRMt, PCRM, PCRMt)
meanFits_long <- xtable(meanFits_long, align = c("l", "l", rep("c", ncol(meanFits_long)-1)),
                        caption="\\raggedright Mean and standard deviation of parameter fits for experiment 2")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(meanFits_long)-1)))
addtorow$command <- c('\\toprule  & \\multicolumn{8}{c}{Model}\\\\  \\cmidrule{2-9} ','[0.2mm]') #'[1mm]'


print(meanFits_long, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(meanFits_long)), booktabs = TRUE,
      caption.placement="top")

print(meanFits_long, type="latex", file="../../Supplement/figures/tablepars2.tex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(meanFits_long)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
#print.xtable()

### Table with BIC comparisons for Supplement
print_BIC_table <- descr_BIC %>% rename(model=model.x) %>%
  mutate(BIC = paste0("$",round(MeanBIC), " (",round(SDBIC), ")$")) %>% 
  select(-MeanBIC, -SDBIC) %>% 
  left_join(select(filter(ungroup(criteria_comparison), criteria=="BIC"), -criteria, -M_diff)) %>%
  mutate(Lower = format(round(Lower, digits=2),nsmall=2), 
         Upper=format(round(Upper, digits=2),nsmall=2),
         BF_power = floor(log(BF_model_vs_WEVmu, base=10)),
         BF = paste0("$", format(round(BF_model_vs_WEVmu/10^BF_power, digits = 1),nsmall=1),
                     "\\times10^{",BF_power, "}$"),
         HDI = paste0("$[", Lower, ", ", Upper, "]$"),
         model= factor(model, levels=c("2DSD", "DDMConf", "dynVib", "IRM", "IRMt", "PCRM", "PCRMt"),
                       labels=c("2DSD",  "DDConf", "dynVis", "IRM", "IRMt", "PCRM", "PCRMt")))%>%
  # labels=c("two-stage dynamical\nsignal detection",
  #          "dynamical visibility",
  #          "independent race model\n(time-dependent confidence)","independent race model\n(Balance of Evidence)",
  #          "drift diffusion confidence",
  #          "anti-correlated race model\n(time-dependent confidence)","anti-correlated race model\n(Balance of Evidence)"))) %>%
  select(-BF_model_vs_WEVmu, -BF_power,-Lower, -Upper)
names(print_BIC_table) <- c("Model", "$M_\\Delta~(SD_\\Delta)$","$BF_{10}$", "95\\%~CI")
print_BIC_table
table_BIC_comparison <- xtable(print_BIC_table, align = c("l", "l", "c", "c", "c"),
                               caption="\\raggedright Results from quantitative mode comparison with dynWEV using BIC for experiment 2")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(table_BIC_comparison)-1)))
addtorow$command <- c('\\toprule  ','[1mm]') #& \\multicolumn{6}{c}{Model}\\\\  \\cmidrule{2-7} 


print(table_BIC_comparison, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BIC_comparison)), booktabs = TRUE,
      caption.placement="top")

print(table_BIC_comparison, type="latex", file="../../Supplement/figures/tableBICs2.tex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BIC_comparison)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
#print.xtable()





### Create a plot for RT-quantiles dependent on choice probability (i.e. discriminability)   ### 

Data_RTQuants_corr_cond <- Data_RTQuants_corr_cond %>%
  mutate(condition = as.numeric(condition),
         X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))
Preds_RTQuants_corr_cond <- Preds_RTQuants_corr_cond %>%
  mutate(condition = as.numeric(condition),
         X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))

Preds_RTQuants_corr_cond<- ungroup(Preds_RTQuants_corr_cond)
Data_RTQuants_corr_cond <- ungroup(Data_RTQuants_corr_cond)
Preds_RTQuants_corr_cond$model1 <- factor(Preds_RTQuants_corr_cond$model, 
                                          levels=c("dynWEV", "2DSD", "dynVib", "DDMConf", "IRMt", "IRM", "PCRMt", "PCRM"),
                                          labels=c("dynamical weighted\nevidence and\nvisibility (dynWEV)",
                                                   "two-stage dynamical\nsignal detection\n(2DSD)",
                                                   "dynamical visibility\n(dynVis)",
                                                   "drift diffusion\nconfidence (DDConf)",
                                                   "independent race\nmodel (time-dependent\nconfidence, IRMt)",
                                                   "independent race\nmodel (Balance of\nEvidence, IRM)",
                                                   "anti-correlated race\nmodel (time-dependent\nconfidence, PCRMt)",
                                                   "anti-correlated race\nmodel (Balance of\nEvidence, PCRM)"))

ggplot()+
  geom_line(data=mutate(Data_RTQuants_corr_cond,, correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))),
            aes(x=p_correct, y=log(q), group=as.factor(p), linetype=correct),color="gray", size=0.7)+  
  geom_line(data=mutate(Preds_RTQuants_corr_cond, correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct")),
  ),#condition = as.factor(condition)
  aes(x=p_correct, y=log(q), group=as.factor(p),color=correct), size=0.7)+
  geom_point( data=mutate(Preds_RTQuants_corr_cond, correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))),
              aes(x=p_correct, y=log(q), color=correct), shape=1,
              size=1.5, fill="white") + 
  geom_point(data=mutate(Data_RTQuants_corr_cond, correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct")),
  ),
  aes(x=p_correct, y=log(q), shape=correct),
  size=1.2, fill="white")+
  scale_color_manual(values= rev(two_colors_correct), breaks=rev(c("Correct", "Wrong")),
                     name = "Predicted",
                     labels=rev(c("Correct", "Wrong"))) +
  scale_shape_manual(values=c(17,21),breaks=rev(c("Correct", "Wrong")),
                     name = "Observed",
                     labels=rev(c("Correct", "Wrong")))  +
  scale_linetype_manual(values=c("dashed", "dashed"),breaks=rev(c("Correct", "Wrong")),
                        name = "Observed",
                        labels=rev(c("Correct", "Wrong")))  +
  scale_x_continuous(name="Choice probability")+
  scale_y_continuous(breaks = log(c(1.5, 2.5, 4, 7)),
                     # labels = paste("log(", c(1.5, 2.5, 4, 7), ")", sep=""), 
                     labels = paste(c(1.5, 2.5, 4, 7), sep=""), 
                     name="Reaction time quantiles [s] (log scaled)")+
  facet_grid(cols=vars(correct), rows=vars(model1), scales = "free_x",
             labeller = )+ #,dir="v"
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        legend.text = element_text(size=9),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(size=7, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9, family = "Times", angle=0),
        strip.text.y.right = element_text(angle = 0),
        legend.box = "horizontal",
        legend.position = "bottom",# legend.position = c(.70, .005), legend.justification = c(0,0),
        legend.direction = "horizontal", legend.spacing.y = unit(0, "lines"),
        legend.margin =margin(0,0,0,1, "cm"), legend.box.spacing = unit(0.2,"lines"),
        #legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.width=unit(1.5,"line"),
        panel.spacing=unit(0, "lines"))

ggsave("figures/RTQuantsCond2.eps",
       width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)

