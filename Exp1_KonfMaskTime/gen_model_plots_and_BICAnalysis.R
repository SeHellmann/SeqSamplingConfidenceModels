####################################################################
###  Generate figures for comparison between model predictions  #### 
###   and empirical data and conduct Bayesian tests for BIC     ####  
####################################################################
#### Structure:                                                 ####
#### 1) Load data and model fits from previous analyses         ####
#### 2) Produce Figures for model fits                          ####
#### 3) BIC Analyses: Bayes t test + Figure for BIC differences ####
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
#### 1) Load data and model fits from previous analyses         ####
load("collected_fitsNpredicts.RData")
levels(Preds_MRating_corr_cond$model) <- c("dynamical weighted evidence\nand visibility", "two-stage dynamical\nsignal detection",
                                           "independent race model\n(time-dependent confidence)","independent race model\n(Balance of Evidence)",
                                           "anti-correlated race model\n(time-dependent confidence)","anti-correlated race model\n(Balance of Evidence)")

### For discussion: Proportion of post-decisional accumulation in response time
Data %>% left_join(select(filter(fits_WEVmodels, model=="WEVmu"), participant, tau)) %>%
  ungroup() %>%
  summarise(prop_post = mean(tau/rt))
Data %>% mutate(easy_error = if_else(as.numeric(condition)>=4 & correct==0, 1, 0)) %>%
  ungroup() %>% summarise(mean(easy_error))


#### 2) Produce Figures for model fits                          ####
####### Figure 5: Plot of Mean Ratings                          ####
pd <- position_dodge(0.02)
Data_MRating_corr_cond <- Data %>% group_by(participant, condition, correct) %>%
  summarise(MRating = mean(rating)) %>%
  summarySEwithin(measurevar = "MRating", idvar = "participant", withinvars = c("condition", "correct")) %>%
  rename(sewithin = se) %>% select(-MRating) %>%
  right_join(mutate(Data_MRating_corr_cond, correct=as.factor(correct)), by=c("condition", "correct"))
  
p_MRating <- ggplot(mutate(Data_MRating_corr_cond, condition = as.integer(condition)),
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  geom_line(data=Preds_MRating_corr_cond, aes(color=as.factor(correct)), size=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.2,position =pd, size=0.4) +
  geom_point(position = pd, fill="white")+
  facet_wrap(.~model, nrow = 3)+ #, dir="v"
  ylab("Mean Confidence")+
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]")+  
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
        legend.position = "bottom", legend.direction = "vertical",#c(.005, .995), legend.justification = c(0,1),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
p_MRating
ggsave("figures/meanRating1.eps", 
       width = 17.62, height=14, units="cm",dpi=1200, device = cairo_ps)


####### Suppl Figure 1: Plot of Ratings Dist.                   ####
ann_text <- expand.grid(correct=c(0,1),
                        condition= levels(Preds_RatingDist_corr_cond$condition ),
                        label=NA, rating=0, p=1)
ann_text$label <- c("A", rep("", 9))
p_ratingdist <- ggplot(data=Data_RatingDist_corr_cond, aes(x=rating, y=p))+
  geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
  geom_point(data=Preds_RatingDist_corr_cond,aes(shape=model, fill=model), position=position_dodge(1), col="black")+
  facet_grid(cols=vars(correct), rows = vars(condition), 
             labeller = labeller(correct=c("0"="Wrong", "1"="Correct"), 
                                 condition=function(x) paste(x, "ms")))+
  scale_x_continuous(
    name = "Confidence [% scale width]", 
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
ggsave("figures/RespDist1.eps", 
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)



####### Figure 6: RTQuantiles accross correct X rating          ####
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
              size=1.5, fill="white") +  
  geom_line(data=Data_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Observed", color="Observed"), size=0.7)+
  geom_line(data=Preds_RTQuants_corr_rating, aes(x=X, y=q, group=as.factor(p), linetype="Model fit", color="Model fit"), size=0.7)+

  geom_point(data=Data_RTQuants_corr_rating, aes(x=X, y=q, shape="Observed"),
             size=1.2)+
  scale_shape_manual(name="class", breaks = c("Model fit", "Observed"), values=c(32, 17))+ 
  scale_x_discrete(name="Confidence x Accuracy", 
                   labels=c("100-80", "60-80","40-60\n Wrong", "20-40",  "0-20", "0-20","20-40", 
                            "40-60\n Correct","60-80", "80-100"))+
  scale_y_continuous(breaks = c(2, 4, 6, 8), name="Reaction Time Quantiles [s]")+
  scale_linetype_manual(name="class", values= c("solid", "dashed"), breaks = c("Model fit", "Observed"))+
  scale_color_manual(name="class", values= two_colors_model, breaks = c("Model fit", "Observed"))+
  facet_wrap(.~model, nrow=3)+#,dir="v"
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
        legend.position = "bottom",# legend.justification = c(1,1),#legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"),
        panel.spacing=unit(0, "lines"))
ggsave("figures/RTQuantsConf1.eps", 
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ungroup(Data) %>% filter(rating==5 & correct==0) %>% summarise(N=n()/nrow(Data))


#### 3) BIC Analyses: Bayes t test + Figure for BIC differences ####

## Compute BIC differences between dynWEV and other models
descr_BIC <- rbind(fits_RMmodels[,c("participant", "model", "BIC")],
                   fits_WEVmodels[,c("participant", "model", "BIC")]) %>%
  filter(model!="WEVmu") %>% 
  left_join(filter(fits_WEVmodels[,c("participant","model", "BIC")], 
                   model=="WEVmu"), by="participant") %>%
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


####### SupplFig 3: Binned BIC diffs btw dynWEV and other models#### 
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
ggsave("figures/BICs1.eps", 
       width=17.62, height=10, dpi=1200, units="cm", device=cairo_ps)
