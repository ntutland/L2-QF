### Look at prescribed fires

library(tidyverse)
library(terra)
library(tidyterra)

getwd()
# setwd("Projects/L2-QF/1.LANDIS-MODEL/LANDIS_run-Fire2-Dry_noFWI/scrapple-fire")
setwd("1.LANDIS-MODEL/LANDIS_run-Fire2-Dry/scrapple-fire")
setwd("..")
FB <- rast("RX_Ignitions_New.tif")
setwd("scrapple-fire")

img_stack <- function(pattern, raster){
  filelist <- list.files(pattern = pattern)
  temp <- rast(filelist[[1]])
  stack <- rast(crs = crs(raster),
              extent = ext(raster),
              resolution = res(raster),
              vals = values(temp))
  names(stack) <- paste0(pattern,"-",1)
  for(i in 2:length(filelist)){
    temp <- rast(filelist[[i]])
    rst <- rast(crs = crs(raster),
                extent = ext(raster),
                resolution = res(raster),
                vals = values(temp))
    names(rst) <- paste0(pattern,"-",i)
    stack <- c(stack,rst)
  }
  return(stack)
}

day_of_fire <- img_stack("day-of-fire",FB)
for(i in 1:50){
  day_of_fire[[names(day_of_fire)[i]]][day_of_fire[[names(day_of_fire)[i]]]>0] <- 1
}
dof_combined <- sum(day_of_fire)
ggplot() + geom_spatraster(data = dof_combined)
writeRaster(dof_combined, "number_of_fires.tif",overwrite=T)
