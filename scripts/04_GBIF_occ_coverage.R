library(spocc)
library(data.table)
library(terra)
library(sf)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(ggpubr)
library(foreach)
sampling_summary <- fread("ClimNiche_database/sampling_summary.csv")

N_raw_occ_gbif <- foreach( sp = sampling_summary$name,.combine = rbind)%do%{
  
  occ_test <- occ( sp ,limit = 1)
  result <- data.table(name = sp,occ_found = occ_test$gbif$meta$found)
  result
  
}

sampling_summary <- merge(sampling_summary,N_raw_occ_gbif[!duplicated(name),c(1,2)])
sampling_summary <- merge(sampling_summary,ClimNicheHub[metric == "mean",])


sampling_summary[,occ_found_log10 := log10(occ_found)]
sampling_summary[,occ_found_log10 := occ_found_log10 - min(occ_found_log10)]

(hist_gbif <- ggplot(sampling_summary,aes(x =occ_found ))+
  geom_histogram(fill = "#76B476",color = "grey20",binwidth = 0.5,boundary = 1,lwd = 1)+
  scale_x_continuous(transform = "log10",
                     breaks = c(100,1000,10000,100000,1000000,10000000),
                     labels = c("100","1,000","10,000","100,000","1,000,000","10,000,000"))+
  theme_classic()+
  labs(x = "Number of Occurrences in GBIF, log10 scale",y = "Number of species"))

ggsave(file.path("figures","GBIF_sampling_figures","hist_GBIF_sample.jpg"),
       hist_gbif,
       dpi = 600,
       unit = "cm",
       width = 18,
       height = 14,
       scale = 0.8)

library(brms)
model_bio10_occ <- brm(bf(bio10 ~ poly(occ_found_log10,2),sigma ~ poly(occ_found_log10,2)),
                       data = sampling_summary,
                       chains = 3,
                       family = gaussian,
                       cores = 3,
                       iter = 3000,
                       backend = "cmdstanr",
                       warmup = 200)
pp_check(model_bio10_occ)
summary(model_bio10_occ)

pred_model <- fitted(model_bio10_occ,sampling_summary,re_formula = NA)
pred_model <- cbind(sampling_summary,pred_model[,1])
pred_model[,Estimate_mean:= V2];pred_model[,V2 := NULL]
pred_model <- cbind(pred_model, fitted(model_bio10_occ,sampling_summary,re_formula = NA,dpar = "sigma")[,1] )
pred_model[,error:= V2];pred_model[,V2 := NULL]
pred_model[,Q2.5:= Estimate_mean - 1.96*error]
pred_model[,Q97.5:= Estimate_mean + 1.96*error]


(model_sampling_plot <- ggplot(sampling_summary,aes(x = occ_found, y = bio10))+
  geom_point(color = "grey20",alpha = 0.75,pch = 16)+
  geom_ribbon(aes( ymin =  Q2.5, ymax = Q97.5),data = pred_model,fill = "#76B476",alpha = 0.25)+
  geom_line(aes(y = Estimate_mean),data = pred_model,color = "#76B476",lwd = 1.5)+
  scale_x_continuous(transform = "log10",
                     breaks = c(1000,10000,100000,1000000),
                     labels = c("1,000","10,000","100,000","1,000,000"))+
  theme_classic()+
  labs( y = "TWQ thermal optimum (°C)", x = "Number of Occurrences in GBIF, log10 scale"))

ggsave(file.path("figures","GBIF_sampling_figures","Model_bio10_sample.jpg"),
       model_sampling_plot,
       dpi = 600,
       unit = "cm",
       width = 18,
       height = 12,
       scale = 0.8)
