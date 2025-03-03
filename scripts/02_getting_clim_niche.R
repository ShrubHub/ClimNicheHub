#### computing climatic niche ####

# This script uses the filtered occurrence and overlay them on the rasters of climatologies
# obtained in CHELSA V2.1. 

#### packages ####

library(data.table)
library(terra)
library(sf)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(ggpubr)

#### meta data and climate data ####

table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"))

all_chelsa_path <- list.files("chelsa_data_V2.1",full.names = T,pattern = ".tif")

chelsa_full <- rast(all_chelsa_path)

names(chelsa_full) <- str_remove(names(chelsa_full),"_1981_2010_V.2.1")
names(chelsa_full) <- str_remove(names(chelsa_full),"CHELSA_")

cell_id <- chelsa_full$bio1
names(cell_id) <- "cell_id"
values(cell_id) <- 1:ncell(cell_id)

chelsa_full <- c(chelsa_full,cell_id)

#### loading occurrences ####

result_occ_path <-  list.files("filtered_occ",full.names = T,pattern = ".RData")

list_of_occ <- lapply(result_occ_path,readRDS)

occ_dt <- rbindlist(list_of_occ);rm(list_of_occ);gc()

occ_dt <- occ_dt[,c(1:6)]

final_count <- occ_dt[,.N,by= name]

enough_occ <- final_count[N>100 ,]

occ_dt <- occ_dt[name%in%enough_occ$name,]

occ_dt_sf <- st_as_sf(occ_dt,coords=c("longitude","latitude"),crs=st_crs(4326))


table_species_name_count <- merge(table_species_name,final_count,
                                  by.x = "ITEX_name",by.y ="name",all.x = T)
hist(final_count$N,nc = 20)
View(final_count)


bugged_sp <- table_species_name_count[is.na(N) & n_occ!= 0,ITEX_name]

#### extraction of climate data ####

all_clim <- extract(chelsa_full,vect(occ_dt_sf))
all_clim <- data.table(all_clim)
occ_dt_clim <- cbind(occ_dt,all_clim[,-"ID"])

#### summary statistics computation ####

col_to_compute <- c(names(chelsa_full)[-33],"longitude","latitude")
col_to_compute <- col_to_compute[c(1,12:19,2:11,20:30,32:34)]  ### ordering, not interested in kg5

round_digit <- c(2,1,3,1,2,2,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,2,0,3,3)

occ_dt_clim[ , metric := "mean" ]
mean_clim_dt <- occ_dt_clim[ , lapply(.SD,mean,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]
mean_clim_dt <- mean_clim_dt[ , mapply(round,.SD,round_digit,SIMPLIFY = F), by = .(name,metric),.SDcols = col_to_compute] ##rounding to the decimal precision of CHELSA

occ_dt_clim[ , metric := "median" ]
med_clim_dt <- occ_dt_clim[ , lapply(.SD,median,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]

occ_dt_clim[ , metric := "Q" ]
Q_clim_dt <- occ_dt_clim[ , lapply(.SD,quantile,probs=c(0.05,0.95),na.rm=T), by = .(name,metric),.SDcols = col_to_compute]
Q_clim_dt[,metric := rep(c("Q05","Q95"),times= nrow(Q_clim_dt)/2)]

Q_clim_dt_tmp <- Q_clim_dt
Q_clim_dt_tmp$metric <- "range"


range_clim_dt <- Q_clim_dt_tmp[,lapply(.SD,diff),by = .(name,metric),.SDcols = col_to_compute];rm(Q_clim_dt_tmp)

get_density_max <- function(x,na.rm=T){
  tmp <- try(density(x,na.rm=na.rm,n = 512*2),silent = T)
  
  if(class(tmp)=="density") return(tmp$x[which.max(tmp$y)]) else return(as.numeric(NA))
  
}
occ_dt_clim[ , metric := "optimum" ]
opt_clim_dt <- occ_dt_clim[ , lapply(.SD,get_density_max,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]
opt_clim_dt <- opt_clim_dt[ , mapply(round,.SD,round_digit,SIMPLIFY = F), by = .(name,metric),.SDcols = col_to_compute] ##rounding to the decimal precision of CHELSA


ClimNicheHub <- rbind(mean_clim_dt,med_clim_dt,Q_clim_dt,range_clim_dt,opt_clim_dt)
ClimNicheHub <- ClimNicheHub[order(name),]

ClimNicheHub[,print(summary(.SD)),by = metric]

#### Export ####

write.table(ClimNicheHub,file.path("ClimNiche_database","climate_summary.csv"),row.names = F,sep = ",")
write.table(table_species_name_count,file.path("ClimNiche_database","sampling_summary.csv"),row.names = F,sep = ",")

table_species_name_count


#### PCA-based approach ####

grid_clim <- occ_dt_clim[!duplicated(cell_id),-c("name","prov","date","key","metric")]
grid_clim <- grid_clim[complete.cases(grid_clim),]
set.seed(0)
sample_of_grid <- sample(grid_clim$cell_id,20000)

plot(grid_clim[sample_of_grid,.(bio1,bio4)])

saveRDS(grid_clim,file.path("Complete_sampling","all_CHELSA_cells.RData"))

saveRDS(occ_dt_clim[,c("name","cell_id")],file.path("Complete_sampling","species_cells.RData"))

occ_dt_clim[,c("name","cell_id")][,.N,by = name]

library(ade4)

to_pca <- grid_clim[cell_id%in%sample_of_grid,-c("longitude" , "latitude","kg5")]
rownames(to_pca) <- to_pca$cell_id
to_pca[,cell_id:=NULL]

to_pca <- to_pca[,-c("gddlgd0","gdgfgd0","fgd","gsp")]

pca_clim_grid <- dudi.pca(to_pca,nf = 4, scannf = F)
inertia.dudi(pca_clim_grid)
s.arrow(pca_clim_grid$co)

col_to_get <- colnames(to_pca)

test_salix <- suprow(pca_clim_grid,occ_dt_clim[name == "Betula pubescens",..col_to_get])

projected_species <- suprow(pca_clim_grid,occ_dt_clim[,..col_to_get])
projected_species <- data.table(projected_species$lisup)
projected_species[,name:=occ_dt_clim$name]

pca_species <- projected_species[,.(Axis1 = mean(Axis1,na.rm = T),Axis2 = mean(Axis2,na.rm = T),
                                    Axis3 = mean(Axis3,na.rm = T),Axis4 = mean(Axis4,na.rm = T)),by=name]

pca_ClimNicheHub <- pca_species

pca_species <- merge(pca_species,ClimNicheHub[metric == "median",])

ggplot(pca_species,aes(x = Axis3, y = Axis4, color = bio10,size = bio4))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()

plot(test_salix$lisup)

write.table(pca_ClimNicheHub,file.path("ClimNiche_database","pca_summary.csv"),row.names = F,sep = ",")

#### misc ####
bootstrap <- foreach()

one_boot <- occ_dt_clim[name == "Achillea millefolium",][sample(1:6904,1000),][ , lapply(.SD,median,na.rm=T), by = name,.SDcols = names(chelsa_full)]

one_boot <- occ_dt_clim[name == "Anthyllis vulneraria",][sample(1:6569,6569,F),]

ggplot(one_boot,aes(x = bio10))+
  theme_bw()+
  geom_histogram(fill="orange",alpha=0.25,color="grey80",mapping = aes(y = after_stat(density)))+
  geom_density(fill="orange",alpha = 0.5)

hist(one_boot$bio1)

occ_dt_2[,lapply(.SD, mean_all),by = name ]

(cols) := lapply(.SD, "*", -1)

density(occ_dt_clim[name == "Salix arctica",bio10],n = 512*2)$x

get_density_max(occ_dt_clim[name == "Salix arctica",bio10])


ggplot(ClimNicheHub[metric == "median",],aes(x = bio10,y = bio12))+
  theme_bw()+
  labs(x = " TWQ °C",y= 'PP mm')+
  geom_point()+
  geom_smooth(method = "lm")

#### Plotting ####
country_shape <- ne_download(50)

plot_one_sp <- function(name_sp, what = "bio10",labs = NULL){
  set.seed(0)
  if(is.null(labs))labs <- what
  
  occ_dt_clim_one_sp <- occ_dt_clim[name ==  name_sp,]
  occ_dt_clim_one_sp_subset <- occ_dt_clim_one_sp
  if(nrow(occ_dt_clim_one_sp)>1000) occ_dt_clim_one_sp_subset <- occ_dt_clim_one_sp[sample(1:nrow(occ_dt_clim_one_sp),1000,F),]
  
  sf_occ <- st_as_sf(occ_dt_clim_one_sp_subset,coords=c("longitude" , "latitude"),crs=st_crs(4326))
  
  bbox <- st_bbox(sf_occ)
  
  print(bbox)
  plot_1 <- ggplot()+
    geom_sf(data = country_shape)+
    geom_sf(data = sf_occ,aes(fill= get(what)),shape=21,color="grey10",alpha = 0.95)+
    scale_fill_gradient(low = "white", high = "orange")+
    coord_sf(xlim = bbox[c(1,3)],ylim = bbox[c(2,4)])+
    theme_bw()+
    labs(fill = labs,title = name_sp)+
    theme(legend.position = "bottom",legend.key.width = unit( 20,"mm"))
 
  plot_2 <- ggplot(occ_dt_clim_one_sp,aes(x = get(what)))+
    theme_bw()+
    labs(x = labs,y="")+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
    geom_density(fill="orange",alpha = 0.5)
  
  final_plot <- ggarrange(plot_1,plot_2,
                          ncol=1,align = "v",heights = c(1.3,1),
                          labels = c(" "," "))
  
  return(final_plot)
     
}

save_a_plot <- function(name_sp, what = "bio10",labs = NULL){
  
  tmp <- plot_one_sp(name_sp,what,labs)
  
  name_sp_2 <- str_replace_all(name_sp," ","_")
  
  file_name <- paste0(name_sp_2,"_",what,".jpg")
    
  ggsave(file.path("figures","individual_species_map",file_name),
         plot = tmp,
         width = 180,
         height = 140,
         unit = "mm",
         dpi = 150)
  
}

save_a_plot("Eriophorum vaginatum","bio10","TWQ (°C)")
save_a_plot("Salix chamissonis","bio10","TWQ (°C)")
save_a_plot("Salix arctica","bio10","TWQ (°C)")
save_a_plot("Anemone narcissiflora","bio10","TWQ (°C)")
save_a_plot("Epilobium palustre","bio10","TWQ (°C)")
save_a_plot("Achillea atrata","bio10","TWQ (°C)")
save_a_plot("Betula pubescens","bio10","TWQ (°C)")
save_a_plot("Betula nana","bio10","TWQ (°C)")
save_a_plot("Androsace obtusifolia","bio10","TWQ (°C)")
save_a_plot("Androsace chamaejasme","bio10","TWQ (°C)")

Androsace obtusifolia
Androsace chamaejasme

Epilobium palustre

Achillea atrata
spocc("")