library(terra)
library(here)
bp <- vect(here("1.LANDIS-MODEL","Shapefiles","burn_plot.shp"))
# writeVector(bp, here("1.LANDIS-MODEL","Shapefiles","burn_plot.kml"))

dem <- rast(here("Hoke_Ground_3ft.tif"))
bp <- project(bp, dem)

new_ext <- ext(bp)
new_ext[1] <- new_ext[1]-((120*3.125))
new_ext[2] <- new_ext[2]+((120*3.125))
new_ext[3] <- new_ext[3]-((120*3.125))
new_ext[4] <- new_ext[4]+((120*3.125))

dem <- crop(dem, new_ext)
dem <- project(dem, "EPSG:32617", res = c(2,2), method = "bilinear")

bp <- vect(here("1.LANDIS-MODEL","Shapefiles","burn_plot.shp"))
new_ext <- ext(bp)
new_ext[1] <- new_ext[1]-100
new_ext[2] <- new_ext[2]+100
new_ext[3] <- new_ext[3]-100
new_ext[4] <- new_ext[4]+100

ext_vect <- vect(new_ext, crs = crs(topo))
plot(topo)
plot(ext_vect, add=T)

topo <- crop(dem, new_ext)

topo <- focal(topo, fun=mean, na.policy="only")
writeRaster(topo, here("Output_Rasters","topo.tif"), overwrite=T)
