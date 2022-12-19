### Code for Supplementary Figure 1 for  Metacogn Sens and fitted tau
library(Hmisc)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(grid)

## Include font style and define colors for plots
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#1f78b4", "#a6cee3", "#33a02c","#b2df8a", "#e31a1c","#fb9a99", "#000000", "#000000")
two_colors_model <- c("#1f78b4", "#b2df8a")

dir.create("figures", showWarnings = FALSE)



load("collected_fitsNpredicts_KonfMaskTime_review.RData")
Data_AUC <- Data %>% group_by(participant) %>%
  summarise(emp_auc= rcorr.cens(rating, correct)[['C Index']],
            .groups = "drop")
metasens_plotdf1 <- Data_AUC %>% 
  left_join(mutate(select(fits_WEVmodels, model, participant, tau), 
                   model=factor(model, levels=c("2DSD", "WEVmu"), 
                                labels=c("2DSD", "dynWEV")))) %>%
  mutate(experiment="Experiment 1")
load("collected_fitsNpredicts_RDK_review.RData")

Data_AUC <- Data %>% group_by(participant) %>%
  summarise(emp_auc= rcorr.cens(rating, correct)[['C Index']],
            .groups = "drop")
metasens_plotdf2 <- Data_AUC %>% 
  left_join(mutate(select(fits_WEVmodels, model, participant, tau), 
                   model=factor(model, levels=c("2DSD", "WEVmu"), 
                                labels=c("2DSD", "dynWEV")))) %>%
  mutate(experiment = "Experiment 2")
ggplot(rbind(metasens_plotdf1, metasens_plotdf2),
       aes(x=emp_auc, y=tau)) +
  geom_point(shape=19, size=3, alpha=0.6)+
  geom_smooth(method="lm", alpha=0)+
  # geom_text(data=Obs_corr, x=.65, y=1.95, hjust=0, vjust=1, size=9, family="Times",
  #           aes(label=rho), col="blue", parse = TRUE)+
  facet_grid(cols=vars(model), rows=vars(experiment), scales = "free")+
  xlab("Observed Type2-AUC")+ylab(expression(Estimated~postdecisional~accumulation~period~tau))+
  theme_bw() +
  theme(text = element_text(size=9, family="Times" ),
        strip.text = element_text(size=9, family = "Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text( angle=90),# switch off minor gridlines
        legend.text.align = 0, legend.text = element_text(size=9),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),#  plot.margin = unit(c(.5,.5,.5,.5), "lines"),
        panel.spacing = unit(.5, "lines"),
        legend.position = c(.825, .25),
        legend.title = element_blank(),
        legend.key.width=unit(1,"line"),
        legend.key.height = unit(1, "line"))
ggsave("figures/AUCTau12.eps",
       width=14, height=10, dpi=1200, units="cm", device=cairo_ps)
Obs_corr <- rbind(metasens_plotdf1, metasens_plotdf2) %>% 
  group_by(model, experiment) %>%
  summarise(rho=paste0("rho","==",round(cor(tau, emp_auc),2)))
Obs_corr
