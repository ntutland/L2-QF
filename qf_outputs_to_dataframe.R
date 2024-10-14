library(here)
library(terra)
library(tidyverse)
library(tidyterra)

#### DEFINE FUNCTIONS ####

output_to_rst <- function(output, out_arr, plot_bounds){
  out_mat <- as.matrix(out_arr)
  out_rst <- rast(nrow = dim(out_mat)[1],
                  ncol = dim(out_mat)[2],
                  crs = crs(plot_bounds),
                  extent = ext(plot_bounds),
                  vals = out_mat,
                  names = output)
  out_rst <- flip(out_rst, direction="vertical")
  return(out_rst)
}

#### ASSEMBLE DATASET ####
runs <- c("Fire2-Dry","Fire5-Dry","Fire2-Wet","Fire5-Wet")
outputs <- c("mass_burnt_pct",
             "surface_consumption",
             "surface_consumption_pct",
             "canopy_consumption",
             "canopy_consumption_pct",
             "max_power",
             "residence_time_power",
             "residence_time_consumption",
             "max_reaction_rate")

burn_domain <- vect(here("1.LANDIS-MODEL",
                         "Shapefiles",
                         "burn_domain.shp"))
burn_plot <- vect(here("1.LANDIS-MODEL",
                         "Shapefiles",
                         "burn_plot.shp"))

for(run in runs){
  run_df = data.frame("Run" = rep(run,89700))
  for(output in outputs){
    out_arr <- read.table(here("7.QUICFIRE-MODEL",
                               "projects", 
                               run,
                               "Arrays",
                               paste0(output,".txt")))
    out_rst <- output_to_rst(output, out_arr, burn_domain)
    writeRaster(out_rst, here("Output_Rasters",paste0(run,"_",output,".tif")), overwrite=T)
    out_rst <- crop(out_rst, burn_plot)
    out_list <- as.vector(values(out_rst))
    out_df <- data.frame(output = out_list)
    run_df <- cbind(run_df, out_df)
    names(run_df)[length(names(run_df))] <- output
  }
  assign(paste0(run,"_df"),run_df)
}



all_data <- bind_rows(`Fire2-Dry_df`, `Fire2-Wet_df`, `Fire5-Dry_df`, `Fire5-Wet_df`)
write.csv(all_data, here("all_data.csv"), row.names = F)

ccsp_list <- list()
for(i in 1:length(runs)){
  df <- read.csv(here("7.QUICFIRE-MODEL",
                      "projects", 
                      runs[i],
                      "canopy_strata.csv"))
  df$run <- runs[i]
  df$strata <- row.names(df)
  ccsp_list[[i]] <- df
}
ccsp_df <- bind_rows(ccsp_list)
write.csv(ccsp_df, here("canopy_strata.csv"), row.names=F)
