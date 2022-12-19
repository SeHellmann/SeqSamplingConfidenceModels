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

temp <- Data %>% group_by(participant, rating, correct) %>%
  summarise(meanRT = mean(rt)) %>%
  summarySEwithin(measurevar = "meanRT", idvar = "participant", withinvars = c("rating", "correct"))%>%
  mutate(correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct"))) %>%
  select(-meanRT)

descr_Data_corr_rating <- Data %>% 
  group_by(participant, rating, correct) %>%
  summarise(MeanRT = mean(rt)) %>%
  group_by(correct, rating) %>%
  summarise(MeanRT = mean(MeanRT)) %>% 
  ungroup() %>%
  mutate(correct=factor(correct, levels=c(0,1), labels=c("Wrong", "Correct")),
         rating = as.factor(rating)) %>%
  left_join(temp, by=c("rating", "correct")) 




descr_Acc <- ggplot(descr_DataAcc,aes(x=condition, y=Accuracy))+
  geom_errorbar(aes(ymin=Accuracy-se, ymax=Accuracy+se), size=0.4)+
  geom_point(size=2)+geom_line(aes(x=as.numeric(condition), y=Accuracy),size=1)+
  scale_x_discrete(name="Coherence [%]")+
  theme_bw()+
  theme(plot.margin = margin(0, 0.3, 0.3, 0.5, "cm"),
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
  theme(plot.margin = margin(0, 0.3, 0.3, 0.5, "cm"),
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

descr_RT_bycond <- ggplot(data = descr_Data_corr,aes(x=condition, y=MeanRT))+
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
  theme(plot.margin = margin(0, 0.3, 0.3, 0.5, "cm"),
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
  annotate("text", label="C", size=14/.pt, x= -0.5, y=RT_ylimits[1]+1.08*diff(RT_ylimits),
           hjust=1, vjust=1, family="Times", fontface="bold")
descr_RT_bycond


RT_ylimits <- c(min(descr_Data_corr_rating$MeanRT - descr_Data_corr_rating$se),
                max(descr_Data_corr_rating$MeanRT + descr_Data_corr_rating$se))
descr_RT_byrating <- ggplot(data = descr_Data_corr_rating,aes(x=rating, y=MeanRT))+
  # geom_point(data = descr_Data,size=2)+
  # geom_line(data = descr_Data,aes(x=as.numeric(rating), y=MeanRT),size=1.3)+
  # geom_errorbar(data = descr_Data, aes(ymin=MeanRT-RT_SER, ymax=MeanRT+RT_SER))+
  geom_line(data=descr_Data_corr_rating, aes(color=as.factor(correct), x=as.numeric(rating), y=MeanRT),size=1)+
  geom_errorbar(data=descr_Data_corr_rating, aes(ymin=MeanRT-se, ymax=MeanRT+se), width=0.5, size=0.4)+
  geom_point(data=descr_Data_corr_rating, aes(shape=as.factor(correct)), size=2, fill="white")+
  scale_x_discrete(name="Confidence", breaks=1:5,
                   labels=c("0-20","20-40", "40-60","60-80", "80-100"))+
  scale_y_continuous(name="Mean response time [s]")+
  scale_shape_manual(breaks=c("Correct", "Wrong"), values=c(21, 17)) + 
  scale_color_manual(breaks=c("Correct", "Wrong"), values=two_colors)+
  theme_bw()+
  theme(plot.margin = margin(0, 0.3, 0.3, 0.5, "cm"),
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
  annotate("text", label="D", size=14/.pt, x= -0.5, y=RT_ylimits[1]+1.08*diff(RT_ylimits), 
           hjust=1, vjust=1, family="Times", fontface="bold")
descr_RT_byrating




Ps <- ggarrange(descr_Acc,descr_Rating,descr_RT_bycond, descr_RT_byrating, 
                nrow=2,ncol=2, common.legend = TRUE, legend="bottom")
Ps
ggsave(plot = Ps,
       file = "figures/descr2.png",
       width = 17.62/1.5, height=11, units = "cm", dpi=1200)
ggsave(plot = Ps,
       file = "figures/descr2.eps",
       width = 17.62/1.5, height=11, units = "cm", dpi=1200, device = cairo_ps)
