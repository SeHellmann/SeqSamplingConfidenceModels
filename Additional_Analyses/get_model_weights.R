### Code to compute model weights from BIC differences between model fits
### 

library(dplyr)
library(tidyr)
library(xtable)
load("collected_fitsNpredicts_KonfMaskTime_review.RData")
fits_WEVmodels$model <- if_else(fits_WEVmodels$model=="2DSD", 'TDSD', 'dynWEV')

names(fits_WEVmodels)
BICweights <- rbind(fits_RMmodels[,c("participant", "model", "BIC", "AIC")],
                    fits_WEVmodels[,c("participant", "model", "BIC", "AIC")],
                    fit_dynVib[,c("participant", "model", "BIC", "AIC")],
                    fit_DDMConf[,c("participant", "model", "BIC", "AIC")]) %>%
  arrange(participant, model) %>%
  group_by(participant) %>%
  mutate(BICdiff = BIC-min(BIC),
         expBICdiff = exp(-0.5*(BIC-min(BIC))),
         wBIC = round((expBICdiff+1e-300) / sum(expBICdiff+1e-300),20)) %>%
  ungroup() %>% 
  select(-c("BIC", "AIC", "BICdiff", "expBICdiff")) %>%
  pivot_wider(values_from=c("wBIC"), 
              names_from= model) #%>% arrange(round(dynWEV, 1), TDSD) 
agg_BICweights <- BICweights %>% summarise(across(.cols =-c("participant"), 
                                                  .fn = mean))%>% 
  mutate(participant = "Mean")
all(agg_BICweights$IRM==0)                                                  
all(agg_BICweights$PCRM==0)
all(agg_BICweights$IRMt==0)
rbind(BICweights, agg_BICweights) %>% select(-IRM, -IRMt, -PCRM)


AICweights <- rbind(fits_RMmodels[,c("participant", "model", "BIC", "AIC")],
                    fits_WEVmodels[,c("participant", "model", "BIC", "AIC")],
                    fit_dynVib[,c("participant", "model", "BIC", "AIC")],
                    fit_DDMConf[,c("participant", "model", "BIC", "AIC")]) %>%
  arrange(participant, model) %>%
  group_by(participant) %>%
  mutate(AICdiff = AIC-min(AIC),
         expAICdiff = exp(-0.5*(AIC-min(AIC))),
         wAIC = round((expAICdiff+1e-300) / sum(expAICdiff+1e-300),3)) %>%
  ungroup() %>% 
  select(-c("BIC", "AIC","AICdiff", "expAICdiff")) %>%
  pivot_wider(values_from=c("wAIC"), 
              names_from= model) #%>% arrange(round(dynWEV, 1), TDSD)
agg_AICweights <- AICweights %>% summarise(across(.cols =-c("participant"),
                                                  .fn=mean)) %>% 
  mutate(participant = "Mean")
all(agg_AICweights$IRM==0)                                                  
all(agg_AICweights$PCRM==0)
all(agg_AICweights$IRMt==0)
rbind(AICweights, agg_AICweights) %>% select(-IRM, -IRMt, -PCRM)

BICweights %>% mutate(across(-participant, .fns = ~round(.x, 1))) %>% 
  summarise(NPCRMt = sum(PCRMt>0.8), 
            NdynWEV = sum(dynWEV>0.8), 
            N2DSD = sum(TDSD >0.8),
            Nmix = sum(dynWEV<1 & dynWEV>=0.2 & TDSD >= 0.2))
agg_BICweights

print_BICweights <- rbind(BICweights, agg_BICweights)
print_BICweights[,2:9] <- format(as.matrix(round(print_BICweights[,2:9], digits = 2), nsmall=2))

print_BICweights <- select(print_BICweights, Participant=participant, dynWEV, "2DSD"=TDSD, DDConf = DDMConf, 
                           dynVis = dynVib, IRM, IRMt, PCRM, PCRMt)
write.csv(print_BICweights, "BIC_weights1.csv", row.names = F)

Ns <- rep(0, 8)
for (i in 1:nrow(print_BICweights)) {
  for (j in 2:ncol(print_BICweights)) {
    if (as.numeric(print_BICweights[i,j])>0.8) {
      Ns[j-1] <- Ns[j-1]+1
      names(Ns)[j-1] <- names(print_BICweights)[j]
      print_BICweights[i,j] <- paste0("\\textbf{", print_BICweights[i,j], "}")
    }
  }
}
Ns
print_BICweights
table_BICweights <- xtable(print_BICweights, align = c("l", "l", "c","c", "c", "c", "c", "c", "c", "c"),
                               caption="\\raggedright BIC weights for each participant and aggregated in experiment 1")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(table_BICweights)-1)), nrow(table_BICweights)-1)
addtorow$command <- c('\\toprule  & \\multicolumn{6}{c}{Model}\\\\  \\cmidrule{2-9}','[0.2mm]', '\\midrule ') # 


print(table_BICweights, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BICweights)), booktabs = TRUE,
      caption.placement="top")

print(table_BICweights, type="latex", file="../../Supplement/figures/tableBICweights1.tex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BICweights)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
#print.xtable()







load("collected_fitsNpredicts_RDK_review.RData")
fits_WEVmodels$model <- if_else(fits_WEVmodels$model=="2DSD", 'TDSD', 'dynWEV')

names(fits_WEVmodels)
BICweights <- rbind(fits_RMmodels[,c("participant", "model", "BIC", "AIC")],
                    fits_WEVmodels[,c("participant", "model", "BIC", "AIC")],
                    fit_dynVib[,c("participant", "model", "BIC", "AIC")],
                    fit_DDMConf[,c("participant", "model", "BIC", "AIC")]) %>%
  arrange(participant, model) %>%
  group_by(participant) %>%
  mutate(BICdiff = BIC-min(BIC),
         expBICdiff = exp(-0.5*(BIC-min(BIC))),
         wBIC = round((expBICdiff+1e-300) / sum(expBICdiff+1e-300),20)) %>%
  ungroup() %>% 
  select(-c("BIC", "AIC", "BICdiff", "expBICdiff")) %>%
  pivot_wider(values_from=c("wBIC"), 
              names_from= model) %>%
  arrange(-round(dynWEV, 1), 
          -round(TDSD,1),
          -round(PCRMt,1), 
          -round(dynVib,1),
          -round(TDSD,3))
agg_BICweights <- BICweights %>% summarise(across(.cols =-c("participant"), 
                                                  .fn = mean))%>% 
  mutate(participant = "Mean")
# all(agg_BICweights$IRM==0)                                                  
# all(agg_BICweights$PCRM==0)
# all(agg_BICweights$IRMt==0)
# rbind(BICweights, agg_BICweights) %>% select(-IRM, -IRMt, -PCRM)


# BICweights %>% mutate(across(-participant, .fns = ~round(.x, 1))) %>% 
#   summarise(NPCRMt = sum(PCRMt>0.8), 
#             NdynWEV = sum(dynWEV>0.8), 
#             N2DSD = sum(TDSD >0.8),
#             Nmix = sum(dynWEV<1 & dynWEV>=0.2 & TDSD >= 0.2))
# agg_BICweights

print_BICweights <- rbind(BICweights, agg_BICweights)
print_BICweights[,2:9] <- format(as.matrix(round(print_BICweights[,2:9], digits = 2), nsmall=2))

print_BICweights <- select(print_BICweights, Participant=participant, dynWEV, "2DSD"=TDSD, DDConf = DDMConf, 
                           dynVis = dynVib, IRM, IRMt, PCRM, PCRMt)
write.csv(print_BICweights, "BIC_weights2.csv", row.names = F)

for (i in 1:nrow(print_BICweights)) {
  for (j in 2:ncol(print_BICweights)) {
    if (as.numeric(print_BICweights[i,j])>0.8) {
      print_BICweights[i,j] <- paste0("\\textbf{", print_BICweights[i,j], "}")
    }
  }
}
print_BICweights

table_BICweights <- xtable(print_BICweights, align = c("l", "l", "c","c", "c", "c", "c", "c", "c", "c"),
                           caption="\\raggedright BIC weights for each participant and aggregated in experiment 2")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(table_BICweights)-1)), nrow(table_BICweights)-1)
addtorow$command <- c('\\toprule  & \\multicolumn{6}{c}{Model}\\\\  \\cmidrule{2-9}','[0.2mm]', '\\midrule ') # 


print(table_BICweights, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BICweights)), booktabs = TRUE,
      caption.placement="top")

print(table_BICweights, type="latex", file="../../Supplement/figures/tableBICweights2.tex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_BICweights)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")

