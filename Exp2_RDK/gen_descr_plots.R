rm(list = ls())
library(Rmisc)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(grid)
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors <- c("#1b9e77", "#fc8d62")

load("JRRDK/collected_fitsNpredicts.RData")
Data_JRRDK <- ungroup(Data)
Data_JRRDK$participant <- Data_JRRDK$participant + 100

load("KonfRDKtime/collected_fitsNpredicts.RData")
Data <- rbind(select(ungroup(Data), participant, condition, correct, rt, Rating, rating),
              select(Data_JRRDK, participant, condition, correct, rt, Rating, rating))
Data <- Data %>%
  mutate(Rating = (Rating+1)/2*100) %>% ungroup()


dir.create("figures", showWarnings = FALSE)
###   Descriptive Plots of Accuracy and RT summaries     #####
descr_DataAcc <- Data %>% group_by(participant, condition) %>%
  summarise(Accuracy = mean(correct)) %>%
  summarySEwithin(measurevar = "Accuracy", idvar = "participant", withinvars = "condition")

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

descr_Acc <- ggplot(descr_DataAcc,aes(x=condition, y=Accuracy))+
  geom_errorbar(aes(ymin=Accuracy-se, ymax=Accuracy+se), size=0.4)+
  geom_point(size=2)+geom_line(aes(x=as.numeric(condition), y=Accuracy),size=1)+
  scale_x_discrete(name="Coherence [%]")+
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
  scale_x_discrete(name="Coherence [%]")+
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
        legend.position = "bottom",#legend.position = "none",
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
  scale_x_discrete(name="Coherence [%]")+
  scale_y_continuous(name="Mean response time [s]", breaks=c(2, 2.5, 3.0, 3.5))+
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
        legend.position = "none", #legend.position = c(.999, 0.999), legend.justification = c(1,1), #legend.position = "bottom"
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.width=unit(2.5,"line"))+
  coord_cartesian(xlim = c(1,5),ylim=RT_ylimits,  clip = 'off')+
  annotate("text", label="C", size=14/.pt, x= -0.5, y=RT_ylimits[1]+1.02*diff(RT_ylimits),
           hjust=1, vjust=1, family="Times", fontface="bold")
descr_RT

Ps <- ggarrange(descr_Acc,descr_Rating,descr_RT, ncol=3, common.legend = TRUE, legend="bottom")
ggsave(plot = Ps,
       file = "figures/descr2.png",
       width = 17.62, height=6, units = "cm", dpi=1200)
ggsave(plot = Ps,
       file = "figures/descr2.eps",
       width = 17.62, height=6, units = "cm", dpi=1200, device = cairo_ps)







###  Descriptive Behavioral measures    ###
# #### Descriptive plots for different sessions              #####
# AccSessions <- Data %>% select(participant, session, condition, correct) %>% 
#   group_by(participant, session, condition) %>%
#   mutate(nrows=n()) %>%
#   summarise(acc = mean(correct),
#             SER = sd(correct)/sqrt(n()))
# pd <- position_dodge(0.1)
# ggplot(AccSessions, aes(x=condition, y=acc, color=as.factor(session), group=session))+
#   geom_line(position=pd)+
#   geom_errorbar(aes(ymin=acc-SER, ymax=acc+SER), colour="black", width=.1,position =pd) +
#   facet_wrap(~participant)
# 
# AccSessions$session <- as.factor(AccSessions$session)
# AccSessions$participant <- as.factor(AccSessions$participant)
# 
# AccSessions$participant <- as.factor(AccSessions$participant)
# 
# anova(lm(data=AccSessions, acc~session))
# anova(lm(data=AccSessions, acc~participant*condition+session))
# anova(lm(data=AccSessions, acc~participant*condition+session),
#       lm(data=AccSessions, acc~participant*condition))
# 
# BFfull <- anovaBF(acc~condition*participant+session, data=AccSessions)
# BFrestr <- anovaBF(acc~condition*participant, data=AccSessions)
# BFrestr/BFfull
#  
# anovaBF(acc~session, data=AccSessions)
# generalTestBF(acc~condition+participant+session, data=AccSessions)
# 
# AccSessions %>% group_by(session) %>%
#   summarise(acc = mean(acc))
