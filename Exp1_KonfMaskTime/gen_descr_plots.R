#rm(list = ls())
if (!exists("Data")) {
  load("collected_fitsNpredicts.RData")
} 
library(Rmisc)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(grid)
windowsFonts(Times=windowsFont("Times New Roman"))
#two_colors <- c("#0072B2","#D55E00")
two_colors <- c("#1b9e77", "#fc8d62")

dir.create("figures", showWarnings = FALSE)
###   Figure 4: Descriptive Plots of Accuracy, Confidence, and RT summaries     #####
temp <- Data %>% group_by(participant, condition) %>%
  summarise(Accuracy = mean(correct)) %>%
  summarySEwithin(measurevar = "Accuracy", idvar = "participant", withinvars = "condition") %>%
  select(-Accuracy)
descr_Data <- Data %>% 
  group_by(participant, condition) %>%
  summarise(Accuracy = mean(correct),
            MeanRating = mean(Rating), 
            MeanRT = mean(rt)) %>%
  group_by(condition) %>%
  summarise(Accuracy = mean(Accuracy), 
            MeanRating = mean(MeanRating), 
            MeanRT = mean(MeanRT)) %>% 
  ungroup() %>% 
  left_join(temp, by="condition")
temp1 <- Data %>% group_by(participant, condition, correct) %>%
  mutate(correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))) %>%
  summarise(MeanRT = mean(rt), 
            MRating = mean(Rating)) %>%
  summarySEwithin(measurevar = "MRating", idvar = "participant", withinvars = c("condition", "correct")) %>%
  rename(seRating = se) %>% select(-MRating)
temp2 <- Data %>% group_by(participant, condition, correct) %>%
  mutate(correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))) %>%
  summarise(MeanRT = mean(rt), 
            MRating = mean(Rating)) %>%
  summarySEwithin(measurevar = "MeanRT", idvar = "participant", withinvars = c("condition", "correct")) %>%
  rename(seRT = se)%>% select(-MeanRT)
descr_Data_corr <- Data %>% 
  group_by(participant, condition, correct) %>%
  summarise(MeanRating = mean(Rating), 
            MeanRT = mean(rt)) %>%
  group_by(condition, correct) %>%
  summarise(MeanRating = mean(MeanRating), 
            MeanRT = mean(MeanRT)) %>% 
  ungroup() %>%
  mutate(correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))) %>%
  left_join(temp1, by=c("condition", "correct")) %>%
  left_join(temp2, by=c("condition", "correct"))


descr_Acc <- ggplot(descr_Data,aes(x=condition, y=Accuracy))+
  geom_errorbar(aes(ymin=Accuracy-se, ymax=Accuracy+se), size = 0.4)+
  geom_point(size=2)+geom_line(aes(x=as.numeric(condition), y=Accuracy),size=1)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]")+
  theme_bw()+
  theme(plot.margin = margin(0, 0.01, 0, 0.5, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = c(.005, 2), legend.justification = c(0,1), #legend.position = "bottom"
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line")) + 
  annotate("text", label="A", size=14/.pt, x= -0.5, y=1.02, hjust=1, vjust=1, family="Times", fontface="bold")+
  coord_cartesian(clip = 'off', xlim = c(1,5), ylim=c(0.5, 1))
descr_Acc

descr_Rating <- ggplot(data = descr_Data_corr,aes(x=condition, y=MeanRating))+
  # geom_point(data = descr_Data,size=2)+
  # geom_line(data = descr_Data,aes(x=as.numeric(condition), y=MeanRating),size=1.3)+
  # geom_errorbar(data = descr_Data, aes(ymin=MeanRating-Rating_SER, ymax=MeanRating+Rating_SER))+
  geom_line(data=descr_Data_corr, aes(color=as.factor(correct), x=as.numeric(condition), y=MeanRating),size=1)+
  geom_errorbar(data=descr_Data_corr, aes(ymin=MeanRating-seRating, ymax=MeanRating+seRating), width=0.5, size=0.4)+
  geom_point(data=descr_Data_corr, aes(shape=as.factor(correct)), size=2, fill="white")+
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]")+
  scale_y_continuous(name="Mean confidence", 
                     #breaks= -1+2/5*0:5,
                     breaks = seq(0, 100, by=20))+
  scale_shape_manual(breaks=c("Correct", "Wrong"), values=c(21, 17)) + 
  scale_color_manual(breaks=c("Correct", "Wrong"), values=two_colors)+
  theme_bw()+
  theme(plot.margin = margin(0, 0.01, 0, 0.5, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom", #legend.position = "none",
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))+
   coord_cartesian(xlim = c(1,5), ylim=c(0, 100),clip = 'off')+  
   annotate("text", label="B", size=14/.pt, x= -0.5, y=102, hjust=1, vjust=1, family="Times", fontface="bold")
descr_Rating


RT_ylimits <- c(min(descr_Data_corr$MeanRT - descr_Data_corr$seRT),
                max(descr_Data_corr$MeanRT + descr_Data_corr$seRT))
descr_RT <- ggplot(data = descr_Data_corr,aes(x=condition, y=MeanRT))+
  # geom_point(data = descr_Data,size=2)+
  # geom_line(data = descr_Data,aes(x=as.numeric(condition), y=MeanRT),size=1.3)+
  # geom_errorbar(data = descr_Data, aes(ymin=MeanRT-RT_SER, ymax=MeanRT+RT_SER))+
  geom_line(data=descr_Data_corr, aes(color=as.factor(correct), x=as.numeric(condition), y=MeanRT),size=1)+
  geom_errorbar(data=descr_Data_corr, aes(ymin=MeanRT-seRT, ymax=MeanRT+seRT), width=0.5, size=0.4)+
  geom_point(data=descr_Data_corr, aes(shape=as.factor(correct)), size=2, fill="white")+
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]")+
  scale_y_continuous(name="Mean response time [s]", breaks=c(2.4, 2.8, 3.2))+
  scale_shape_manual(breaks=c("Correct", "Wrong"), values=c(21, 17)) + 
  scale_color_manual(breaks=c("Correct", "Wrong"), values=two_colors)+
  theme_bw()+
  theme(plot.margin = margin(0, 0.01, 0, 0.5, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position="none", #legend.position = c(.005, 0.995), legend.justification = c(0,1), #legend.position = "bottom"
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))+
  coord_cartesian(xlim = c(1,5),ylim=RT_ylimits,  clip = 'off')+  
  annotate("text", label="C", size=14/.pt, x= -0.5, y=RT_ylimits[1]+1.02*diff(RT_ylimits), 
           hjust=1, vjust=1, family="Times", fontface="bold")
descr_RT

# Ps <- grid.arrange(descr_Acc,descr_Rating,descr_RT, layout_matrix=rbind(c(1,2,3)))
# #annotate_figure(Ps, fig.lab = c("A","B","C"), fig.lab.pos = "top.left")

Ps <- ggarrange(descr_Acc,descr_Rating,descr_RT, ncol=3, common.legend = TRUE, legend="bottom")
ggsave(plot = Ps,
       file = "figures/descr1.eps",
       width = 17.62, height=6, units = "cm", dpi=1200, device = cairo_ps)
