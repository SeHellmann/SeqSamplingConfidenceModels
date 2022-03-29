rm(list = ls())
library(tidyverse)
#install.packages("xtable")
library(xtable)

load("collected_fitsNpredicts.RData")
levels(Preds_MRating_corr_cond$model) <- c("dynamical weighted evidence\nand visibility", "two-stage dynamical\nsignal detection",
                                           "independent race model\n(time-dependent confidence)","independent race model\n(Balance of Evidence)",
                                           "anti-correlated race model\n(time-dependent confidence)","anti-correlated race model\n(Balance of Evidence)")
Fits <- bind_rows(fits_WEVmodels, rename(fits_RMmodels, A=a, B=b))
meanFits <- Fits %>%
  select(-c("participant", "k", "N", "AICc", "negLogLik")) %>%
  group_by(model) %>%
  summarise(across(.cols = everything(),
                   ~paste(round(mean(.x[abs(.x)<1e+24], na.rm=TRUE),2),
                          " (", round(sd(.x[abs(.x) <1e+24], na.rm = TRUE), 2), ")", sep=""))) %>%
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


### Use this to produce a wide table
# meanFits <- column_to_rownames(meanFits, "model")
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="../Supplement/figures/tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)

### Use this to produce a long table
meanFits_long <- meanFits %>% pivot_longer(cols = 2:ncol(meanFits), names_to = "Parameter") %>%
  pivot_wider(id_cols = "Parameter", names_from = "model", values_from = value) %>%
  select(Parameter, dynWEV = WEVmu, '2DSD', IRM, IRMt, PCRM, PCRMt)
meanFits_long <- xtable(meanFits_long, align = c("l", "l", rep("c", ncol(meanFits_long)-1)),
                        caption="\\raggedright Mean and standard deviation of parameter fits for experiment 1")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(meanFits_long)-1)))
addtorow$command <- c('\\toprule  & \\multicolumn{6}{c}{Model}\\\\  \\cmidrule{2-7} ','[1mm]')


print(meanFits_long, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(meanFits_long)), booktabs = TRUE,
      caption.placement="top")

print(meanFits_long, type="latex", file="../Supplement/figures/tablepars1.tex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(meanFits_long)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
#print.xtable()
