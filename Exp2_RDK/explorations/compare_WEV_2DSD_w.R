library(tidyverse)
library(BayesFactor)
load("../collected_fitsNpredicts_bothRDKExperiments.RData")
X_IRM <- (fits_RMmodels %>% filter(model=="IRMt"))[["BIC"]]
X_PCRM <- (fits_RMmodels %>% filter(model=="PCRMt"))[["BIC"]]
Y_WEV <- (fits_WEVmodels %>% filter(model=="WEVmu"))[["BIC"]]

temp <- fits_WEVmodels %>% mutate(model=if_else(model=="WEVmu", "WEVmu", "TDSD")) %>%
  pivot_wider(id_cols = c("participant"), names_from = model, values_from = BIC) %>%
  left_join(select(filter(fits_WEVmodels, !is.na(w)), c("participant", "w"))) %>%
  mutate('TwoDSD_better_than_WEVmu' = (WEVmu > TDSD),
         BICDiff = WEVmu - TDSD)
ggplot(temp) +
  geom_boxplot(aes(x=TwoDSD_better_than_WEVmu, y=w))
ggplot(temp) +
  geom_point(aes(x=w, y=BICDiff, color=TwoDSD_better_than_WEVmu))+
  theme(legend.position = "bottom")
max((subset(temp, TwoDSD_better_than_WEVmu==FALSE))$w)
min((subset(temp, TwoDSD_better_than_WEVmu!=FALSE))$w)