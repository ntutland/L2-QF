library(tidyverse)
library(here)
library(GGally)
library(terra)
library(ggthemes)
library(tidyterra)

dat <- read.csv(here("all_data.csv"))

dat_summ <- dat %>%
  group_by(Run) %>%
  summarize(total_mass_burnt = sum(mass_burnt_pct)/n(),
            total_surf_cons = sum(surface_consumption, na.rm=T),
            total_surf_cons_pct = sum(surface_consumption_pct, na.rm=T)/n(),
            total_canopy_cons = sum(canopy_consumption, na.rm=T),
            total_canopy_cons_pct = sum(canopy_consumption_pct, na.rm=T)/n(),
            average_max_power = mean(max_power[max_power > 0], na.rm=T),
            average_res_time_power = mean(residence_time_power[residence_time_power > 0], na.rm=T),
            average_res_time_cons = mean(residence_time_consumption[residence_time_consumption > 0], na.rm=T),
            average_max_react_rate = mean(max_reaction_rate[max_reaction_rate > 0], na.rm=T)) %>%
  pivot_longer(cols = 2:10,
               names_to = "variable",
               values_to = "value")

dat_sd <- dat %>%
  group_by(Run) %>%
  summarize(average_max_power = sd(max_power[max_power > 0], na.rm=T),
            average_res_time_power = sd(residence_time_power[residence_time_power > 0], na.rm=T)) %>%
  pivot_longer(cols = 2:3,
               names_to = "variable",
               values_to = "sd")

dat_summ_sd <- left_join(dat_summ, dat_sd, by = join_by(Run,variable))

dat_summ_sd$Rotation <- rep(NA, nrow(dat_summ_sd))
dat_summ_sd$Climate <- rep(NA, nrow(dat_summ_sd))
for(i in 1:nrow(dat_summ_sd)){
  dat_summ_sd$Rotation[i] <- strsplit(dat_summ_sd$Run[i],"-")[[1]][1]
  dat_summ_sd$Climate[i] <- strsplit(dat_summ_sd$Run[i],"-")[[1]][2]
}

dat_summ_clean <- dat_summ_sd %>%
  select(Rotation, Climate, variable, value, sd) %>%
  filter(variable %in% c("total_surf_cons_pct",
                         "total_canopy_cons_pct",
                         "average_max_power",
                         "average_res_time_power")) %>%
  mutate(Rotation = if_else(Rotation== "Fire2", "2-Year", "5-Year"),
         Climate = if_else(Climate=="Dry", "Hot-Dry", "Hot-Wet"))

my_colors = colorblind_pal()(8)[2:8]

dat_cons <- dat_summ_clean %>% 
  filter(variable %in% c("total_surf_cons_pct","total_canopy_cons_pct")) %>%
  select(-sd) %>%
  mutate(variable = case_when(variable=="total_surf_cons_pct" ~ "Total Surface\nConsumption (%)",
                              variable=="total_canopy_cons_pct" ~ "Total Canopy\nConsumption (%)"))

consumption_plot <- dat_cons %>% ggplot() +
  geom_bar(stat="identity",aes(Rotation,value, fill=Climate),position=position_dodge()) +
  facet_wrap(.~variable, scales="free_y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Fire Rotation",
       y = "",
       fill = "Climate\nScenario") +
  theme_bw()

ggsave("consumption_comparison.jpg", consumption_plot, path=here("Plots"),height=8,width=16,units = "cm")

dat_surf_eng <- dat_summ_clean %>%
  filter(variable %in% c("average_max_power","average_res_time_power")) %>%
  mutate(variable = case_when(variable=="average_max_power" ~ "Average Maximum\nPower at Surface (kW/s)",
                              variable=="average_res_time_power" ~ "Average Fire Residence\nTime at Surface (s)")) %>%
  mutate(sd_plus = value+sd,
         sd_minus = value-sd)

dat_surf_eng %>% ggplot(aes(x=Rotation,ymin=sd_minus,ymax=sd_plus,fill=Climate)) +
  geom_bar(aes(y=value),stat="identity",position=position_dodge()) +
  geom_errorbar(position = position_dodge(0.9), width=0.5) +
  facet_wrap(.~variable, scales="free_y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Fire Rotation",
       y = "",
       fill = "Climate\nScenario") +
  theme_bw()


#####################
burn_plot <- vect(here("1.LANDIS-MODEL","Shapefiles","burn_plot.shp"))
burn_domain <- vect(here("1.LANDIS-MODEL","Shapefiles","burn_domain.shp"))

f2d_canopy_cons_pct <- rast(here("Output_Rasters","Fire2-Dry_canopy_consumption_pct.tif"))
f2w_canopy_cons_pct <- rast(here("Output_Rasters","Fire2-Wet_canopy_consumption_pct.tif"))
f5d_canopy_cons_pct <- rast(here("Output_Rasters","Fire5-Dry_canopy_consumption_pct.tif"))
f5w_canopy_cons_pct <- rast(here("Output_Rasters","Fire5-Wet_canopy_consumption_pct.tif"))
canopy_consumption <- c(f2d_canopy_cons_pct,f2w_canopy_cons_pct,f5d_canopy_cons_pct,f5w_canopy_cons_pct)
names(canopy_consumption) <- c("Fire2-Dry","Fire2-Wet","Fire5-Dry","Fire5-Wet")

canopy_consumption <- project(canopy_consumption, crs(burn_domain))
ext(canopy_consumption) <- ext(burn_domain)
canopy_consumption <- crop(canopy_consumption, burn_plot)

canopy_cons_map <- ggplot() +
  geom_spatraster(data=canopy_consumption) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr) +
  labs(fill = "Canopy\nConsumption\n(%)") +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave("canopy_cons_maps.jpg", canopy_cons_map, path=here("Plots"),height=14,width=16,units = "cm")

#
f2d_surface_cons_pct <- rast(here("Output_Rasters","Fire2-Dry_surface_consumption_pct.tif"))
f2w_surface_cons_pct <- rast(here("Output_Rasters","Fire2-Wet_surface_consumption_pct.tif"))
f5d_surface_cons_pct <- rast(here("Output_Rasters","Fire5-Dry_surface_consumption_pct.tif"))
f5w_surface_cons_pct <- rast(here("Output_Rasters","Fire5-Wet_surface_consumption_pct.tif"))
surface_consumption <- c(f2d_surface_cons_pct,f2w_surface_cons_pct,f5d_surface_cons_pct,f5w_surface_cons_pct)
names(surface_consumption) <- c("Fire2-Dry","Fire2-Wet","Fire5-Dry","Fire5-Wet")

surface_consumption <- project(surface_consumption, crs(burn_domain))
ext(surface_consumption) <- ext(burn_domain)
surface_consumption <- crop(surface_consumption, burn_plot)

surface_cons_pct_map <- ggplot() +
  geom_spatraster(data=surface_consumption) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr) +
  labs(fill = "Surface\nConsumption\n(%)") +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave("surface_cons_pct_maps.jpg", surface_cons_pct_map, path=here("Plots"),height=14,width=16,units = "cm")

#
f2d_surface_cons <- rast(here("Output_Rasters","Fire2-Dry_surface_consumption.tif"))
f2w_surface_cons <- rast(here("Output_Rasters","Fire2-Wet_surface_consumption.tif"))
f5d_surface_cons <- rast(here("Output_Rasters","Fire5-Dry_surface_consumption.tif"))
f5w_surface_cons <- rast(here("Output_Rasters","Fire5-Wet_surface_consumption.tif"))
surface_consumption <- c(f2d_surface_cons,f2w_surface_cons,f5d_surface_cons,f5w_surface_cons)
names(surface_consumption) <- c("Fire2-Dry","Fire2-Wet","Fire5-Dry","Fire5-Wet")

surface_consumption <- project(surface_consumption, crs(burn_domain))
ext(surface_consumption) <- ext(burn_domain)
surface_consumption <- crop(surface_consumption, burn_plot)

surface_cons_map <- ggplot() +
  geom_spatraster(data=surface_consumption) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr) +
  labs(fill = "Surface\nConsumption\n(g/m^3)") +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave("surface_cons_maps.jpg", surface_cons_map, path=here("Plots"),height=14,width=16,units = "cm")
