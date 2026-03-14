



ClimNicheHub_test <- fread("ClimNiche_database/ClimNiche_climate_summary.csv")
ClimNicheHub_test <- ClimNicheHub_test[metric == "Q95",]

Clim_Nich_pnaspaper <- fread('STI_plant.csv')

lancaster_tolerance <- fread("lancaster_tolerance.csv")
lancaster_tolerance[,name:=str_replace(Plant_species_clean,"_"," ")]




lancaster_tolerance

ClimNicheHub_test <- merge(ClimNicheHub_test,Clim_Nich_pnaspaper,by.x = "name", by.y = "Species")

ClimNicheHub_test <- merge(ClimNicheHub_test,lancaster_tolerance,by ="name")
ClimNicheHub_test <- ClimNicheHub_test[Tolerance_measure == "LT50" & Tmin_or_Tmax == "Tmax",]

table(ClimNicheHub_test$Tolerance_measure,ClimNicheHub_test$Tmin_or_Tmax)


ggplot(ClimNicheHub_test,aes(x = STI, y = bio1))+
  theme_bw()+
  geom_point()+
  geom_abline(slope = 1,intercept = 0,lty = 2, color = "grey50")+
  geom_smooth(method = "lm")

ggplot(ClimNicheHub_test,aes(x = bio1, y = Thermal_Tolerance_deg.C))+
  theme_bw()+
  geom_point()+
  geom_abline(slope = 1,intercept = 0,lty = 2, color = "grey50")+
  geom_smooth(method = "lm")




lm_niches <- lm(bio1 ~ STI,data = ClimNicheHub_test)
summary(lm_niches)

ClimNicheHub_test[,cor(bio1,STI,method = "spearman")]
