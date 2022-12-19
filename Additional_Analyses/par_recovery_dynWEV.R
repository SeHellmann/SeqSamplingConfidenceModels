##### Script to check parameter recovery for 2DSD and dynWEV

library(dynConfiR)
library(tidyverse)
library(ggpubr)
library(grid)
library(gridExtra)

load("../Exp1_KonfMaskTime/collected_fitsNpredicts.RData")
pars_Exp1 <- fits_WEVmodels %>% filter(model=="WEVmu") %>%
  mutate(experiment="Experiment 1")
Ns_Exp1 <- Data %>% group_by(participant) %>%
  summarise(Ntrials = n())
load("../Exp2_RDK/collected_fitsNpredicts_bothRDKExperiments.RData")
pars_Exp2 <- fits_WEVmodels %>% filter(model=="WEVmu") %>%
  mutate(experiment="Experiment 2",
         participant = participant+ 100)
Ns_Exp2 <- Data %>% group_by(participant) %>%
  summarise(Ntrials = n()) %>%
  mutate(participant = participant + 100)
true_pars <- rbind(pars_Exp1, pars_Exp2)
Ns <- rbind(Ns_Exp1, Ns_Exp2)

if (!file.exists("par_recovery/recovered_parameters.RData")) {
  simulate_model <- function(params) {
    model <- params[1]
    cur_N <- ceiling(filter(Ns, participant==params[["participant"]])[["Ntrials"]]/10)
    parnames <- names(params)[c(3:22, 29:31)]
    params <- as.data.frame(t(as.numeric(params[parnames])))
    names(params) <- parnames
  
    sims <- simulateRTConf(params, n=cur_N, model="dynWEV", simult_conf = TRUE)
  }
  set.seed(992)
  #gathered_simus <- apply(true_pars, 1, simulate_model)
  gathered_simus <- data.frame()
  for (i in 1:nrow(true_pars)) {
    params <- true_pars[i,]
    sims <- simulate_model(params)
    sims$participant <- params$participant
    sims$model <- params$model
    sims$experiment <- params$experiment
    gathered_simus <- rbind(gathered_simus, sims)
  }
  
  dir.create("par_recovery", showWarnings = FALSE)
  save(gathered_simus, true_pars, file="par_recovery/saved_simulations.RData")
  
  
  fit_dynWEV_recovery <- fitRTConfModels(gathered_simus, "dynWEV", nRatings = 5, fixed = list(sym_thetas=FALSE),
                                         restr_tau = "simult_conf", 
                                         logging=TRUE,opts = list(nAttempts=4),
                                         parallel="both", n.cores = c(6, 4))
  
  save(fit_dynWEV_recovery, 
       file="par_recovery/recovered_parameters.RData")

}



load("par_recovery/recovered_parameters.RData")


fit_dynWEV_recovery$class <- "est"
fit_dynWEV_recovery$experiment <- if_else(fit_dynWEV_recovery$participant >= 100, "Experiment 2", "Experiment 1")
true_pars$class <- "true"
true_pars <- true_pars %>% rename(svis=sig, sigvis = sigmu)
collected_par_recovery <- fit_dynWEV_recovery %>% select(-setdiff(names(fit_dynWEV_recovery), names(true_pars))) %>%
  rbind(true_pars) %>%
  pivot_longer(cols = t0:sigvis, 
               names_to = "parameter", values_to = "value") %>%
  arrange(participant, model, parameter, class)
collected_par_recovery[is.na(collected_par_recovery$parameter),]
head(collected_par_recovery)

collected_par_recovery_long <- collected_par_recovery %>% 
  pivot_wider(id_cols = c("participant", "parameter", "experiment"), names_from = class, values_from = value)

par_labels <- c('nu[1]', 'nu[2]', 'nu[3]', 'nu[4]', 'nu[5]', 's[nu]', 'a', 'z', 's[z]', 
                't[0]', 's[t0]', 'tau', 'w', 'sigma[V]', 's[V]',
                paste0('theta[', -1,'~', 1:4,']'),#'',
                paste0('theta[', 1,'~',  1:4,']'))
est_cors <- subset(collected_par_recovery_long, abs(est) < 1e+6 & abs(true) < 1e+6) %>%
  group_by(parameter) %>%
  summarise(cor = cor(est, true),
            x_min=min(true), 
            y_max=max(est)) %>%
  mutate(parameter = factor(parameter, levels = 
                              c('v1', 'v2', 'v3', 'v4', 'v5', 'sv', 'a', 'z', 'sz',
                                't0', 'st0', 'tau', 'w', 'sigvis', 'svis',
                                'thetaLower1', 'thetaLower2', 'thetaLower3', 'thetaLower4',# '',
                                'thetaUpper1', 'thetaUpper2', 'thetaUpper3', 'thetaUpper4'
                              ),labels = par_labels 
                              ))
ggplot(est_cors, aes(x=parameter, y=cor))+ #, shape=model, col=model
  geom_point(size=4)+
  scale_x_discrete(labels=parse(text=levels(est_cors$parameter)))+
  theme_bw()+ylim(c(-0.5, 1))






plot_recovery <- subset(collected_par_recovery_long, abs(est) < 1e+6 & abs(true) < 1e+6) %>%
  left_join(Ns) %>% left_join(est_cors) %>%
  arrange(cor)%>%
  mutate(parameter = factor(parameter, levels = 
                              c('v1', 'v2', 'v3', 'v4', 'v5', 'sv', 'a', 'z', 'sz',
                                't0', 'st0', 'tau', 'w', 'sigvis', 'svis',
                                'thetaLower1', 'thetaLower2', 'thetaLower3', 'thetaLower4', #'',
                                'thetaUpper1', 'thetaUpper2', 'thetaUpper3', 'thetaUpper4'
                              ),labels = par_labels ))

windowsFonts(Times=windowsFont("Times New Roman"))
#two_colors <- c("#0072B2","#D55E00")
# two_colors <- c("#1b9e77", "#fc8d62")

dir.create("figures", showWarnings = FALSE)

p_rec <- ggplot(plot_recovery, aes(x=true, y=est, shape=experiment))+
  geom_smooth(data=plot_recovery, method="lm", alpha=0.5, size=0.9, 
              aes(x=true, y=est), formula = 'y~x', inherit.aes = FALSE)+
  geom_abline(col="black", size=1, alpha=0.2)+
  geom_point(alpha=0.5)+ scale_shape_discrete("")+
  geom_text(data = mutate(est_cors, experiment="Experiment 1"),
            aes(x=x_min, y=y_max,
                label=paste0(format(round(cor, 2), nsmall=2))),#"rho == ",
            color="black",
            parse=TRUE, hjust=0, vjust=1)+
  facet_wrap(.~parameter, scales = "free", labeller = label_parsed, 
             drop = FALSE, ncol=4)+
  theme_bw()+ylab("Recovered parameter value")+xlab("True parameter value")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = c(0.78, 0.15), legend.justification = c(0, 1),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9))
p_rec
ggsave(p_rec, file = "figures/par_recovery.eps",
              width = 17.62, height=22.6, units = "cm", dpi=1200, device = cairo_ps)
ggsave(p_rec, file = "../../Supplement/figures/results/par_recovery.eps",
       width = 17.62, height=22.6, units = "cm", dpi=1200, device = cairo_ps)

# g <- ggplotGrob(p_rec)
# # get the grobs that must be removed
# rm_grobs <- g$layout$name %in% c("panel-4-5","strip-t-5-4")
# # remove grobs
# g$grobs[rm_grobs] <- NULL
# g$layout <- g$layout[!rm_grobs, ]
# ## move axis closer to panel
# g$layout[g$layout$name == "axis-b-5-4", c("t", "b")] = c(14.5, 14.5)
# grid.newpage()
# grid.draw(g)
# 
# ggsave(plot = g,
#        file = "figures/par_recovery.eps",
#        width = 17.62, height=22.6, units = "cm", dpi=1200, device = cairo_ps)
# ggsave(plot = g,
#        file = "figures/par_recovery_low.eps",
#        width = 17.62, height=22.6, units = "cm", dpi=600, device = cairo_ps)




df <- subset(collected_par_recovery_long, abs(est) > 1e+6 | abs(true) > 1e+6) 
cor(fit_dynWEV_recovery$sigvis, fit_dynWEV_recovery$svis)
