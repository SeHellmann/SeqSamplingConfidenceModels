## Code for Figure 3 with simulations of confidence in the dynWEV model
## 
## 
library(dynConfiR)
library(tidyverse)
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#1f78b4", "#a6cee3", "#33a02c","#b2df8a", "#e31a1c","#fb9a99", "#000000", "#000000")
two_colors_model <- c("#1f78b4", "#b2df8a")

set.seed(2201)
paramDf <- expand.grid(v1 = 0, v2=1, v3=2, v4=3, v5=4, v6=5, sv=0.5, 
                       a=1, z=0.5, sz=0, t0=0, st0=0, 
                       sigvis=1, svis=1, tau=1, 
                       w=c(0.1, 0.4, 0.6, 0.9),
                       theta1=0, theta2=0.5, theta3=1, theta4=1.5,
                       model="dynWEV")

#simus <- simulateWEV(paramDf, n=5000, model="dynWEV", simult_conf = TRUE, method = "Cpp")
simus <- paramDf %>% group_by(w) %>%
  summarise(simulateRTConf(cur_data_all(), n=10^5))
meanConf <- simus %>% group_by(w, condition, correct) %>% 
  summarise(meanConf = mean(rating))


unique(simus$stimulus)
nrow(simus) / 6 / 2 / 5

#pd <- position_dodge(0.02)
p_MRating <- ggplot(mutate(meanConf, condition = as.integer(condition)),
                    aes(x=condition, y=meanConf, group = as.factor(correct), color=as.factor(correct))) +
  geom_line(data=meanConf, size=1.5)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  #geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.25,position =pd, size=0.6) +
  geom_point(fill="white", size=1.8)+
  facet_wrap(.~w, nrow=2)+ #dir="v"
  ylab("Mean confidence")+
  scale_x_continuous(name=expression(Mean~drift~rate~nu), 
                     breaks=1:6, labels=0:5)+
  scale_y_continuous(breaks=1:5, limits = c(1,5))+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "",
                     labels=c("Correct", "Wrong")) +
  # scale_shape_manual(values=c(21,17),breaks=c(1,0),
  #                    name = "Observed",
  #                    labels=c("Correct", "Wrong"))  +
  #scale_y_continuous(limits = c(1, 5), breaks=1:5, labels = c("0", "25", "50", "75", "100"), name="Mean Confidence")+
  theme_bw() +
  #guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position="bottom",#legend.position = c(.01, .995), legend.justification = c(0,1),
        legend.direction = "horizontal", #legend.position = "bottom"
        legend.title = element_blank(),
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.margin = margin(0, 0, 0, 0, "cm"),#legend.direction = "horicontal",
        legend.box.margin = margin(-0.3, 0, 0, 0, "cm"),
        legend.key.width=unit(1.5,"line"),
        panel.spacing=unit(0, "lines"))
p_MRating
# ggsave("figures/simul_dynWEV.eps",
#        width = 17.62/2, height=12, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/simul_dynWEV.eps",
       width = 17.62/2, height=8, units="cm",dpi=1200, device = cairo_ps)

