getwd()
# setwd("Projects/L2-QF")
setwd("..")

library(terra)
library(tidyverse)
wkt <- crs(rast("1.LANDIS-MODEL/LANDIS_run-Fire2-Dry/IC_Map_FortBragg.tif"))
extn <- ext(rast("1.LANDIS-MODEL/LANDIS_run-Fire2-Dry/IC_Map_FortBragg.tif"))
aoi <- vect("1.LANDIS-MODEL/Shapefiles/burn_plot.shp")
runlist <- c("Fire2-Dry", "Fire2-Wet", "Fire5-Dry", "Fire5-Wet")
# runlist <- c("Fire2-Dry", "Fire2-Wet")
# runlist <- c("LANDIS_run")

df_list <- list()
for(run in runlist){
  needle_mat <- matrix(ncol = 6, nrow = 50)
  litter_mat <- matrix(ncol = 6, nrow = 50)
  for(i in 1:50){
    litter <- rast(paste0("1.LANDIS-MODEL/","LANDIS_run-",run,"/NECN/SurfaceLitterBiomass-",i,".img"))
    litter[litter==0] <- NA
    crs(litter) <- wkt
    ext(litter) <- extn
    litter <- mask(crop(litter,aoi),aoi)
    l_val <- values(litter)[,1]
    for(j in 1:6){
      l_smm <- as.numeric(summary(l_val)[j])
      litter_mat[i,j] <- l_smm
    }
  }
  litter_df <- as.data.frame(litter_mat)
  litter_df$year <- as.numeric(row.names(litter_df))
  names(litter_df) <- c("min","q1","med","mean","q3","max","year")
  df_list[[run]] <- litter_df
}
surf_df <- bind_rows(df_list, .id = "run")

surf_long <- surf_df %>%
  select(-q1,-med,-q3) %>%
  pivot_longer(cols = c("min","mean","max"),
               names_to = "stat",
               values_to = "value")

surf_plot <- surf_long %>%
  ggplot(aes(year,value, color = stat))+
  geom_point() +
  geom_line() +
  facet_wrap(.~run) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Year",
       y = "Biomass (g/m3)",
       color = "Summary\nStatistic") +
  theme_bw()
surf_plot
ggsave("Surface_Fuels_FtBragg.jpg",surf_plot,height = 5, width = 12,units = "in")
