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
  mutate(Rotation = case_when(Rotation == "Fire2" ~ "2-Year", 
                              Rotation == "Fire5" ~ "5-Year",
                              TRUE ~ "No Rx Fire"),
         Climate = if_else(Climate=="Dry" | is.na(Climate), "Hot-Dry", "Hot-Wet"))

my_colors = colorblind_pal()(8)[2:8]
my_colors_rev = colorblind_pal()(8)[c(3,2)]

dat_cons <- dat_summ_clean %>% 
  filter(variable %in% c("total_surf_cons_pct","total_canopy_cons_pct")) %>%
  select(-sd) %>%
  mutate(variable = case_when(variable=="total_surf_cons_pct" ~ "Surface",
                              variable=="total_canopy_cons_pct" ~ "Canopy"),
         Climate = factor(Climate, levels = c("Hot-Wet","Hot-Dry")))

consumption_bar <- dat_cons %>% 
  mutate(Rotation = factor(Rotation, levels = c("5-Year","2-Year"))) %>%
  ggplot() +
  geom_bar(stat="identity",aes(Rotation,value, fill=Climate),position=position_dodge2(preserve="single")) +
  facet_wrap(.~variable, ncol=1, strip.position = "right") +
  scale_y_continuous(limits = c(0,102), expand=c(0,0)) +
  scale_fill_manual(values = my_colors_rev, guide = guide_legend(reverse=T)) +
  labs(x = "Fire Rotation",
       y = "Total Fine Fuel Consumption (%)",
       fill = "Climate\nScenario") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.8,0.77),
        legend.background = element_rect(color = "black"))

ggsave("consumption_comparison.jpg", consumption_bar, path=here("Plots"),
       height=3,width=6, units = "in")

consumption_point <- dat_cons %>% 
  ggplot() +
  geom_point(aes(Rotation,value, color=Climate, shape=variable)) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Fire Rotation",
       y = "Fine Fuel Consumption (%)",
       color = "Climate\nScenario",
       shape = "Fuel\nType") +
  theme_bw()


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
  scale_fill_gradient(low = "lightgreen", high = "black", na.value = "white") +
  facet_wrap(~lyr, ncol =4) +
  labs(fill = "Canopy\nConsumption\n(%)") +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "bottom")
canopy_cons_map

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

###################
# Canopy consumption by strata

ccsp <- read.csv(here("canopy_strata.csv"))
ccsp %>%
  ggplot() +
  geom_line(aes(strata, canopy_consumption_strata, color=run)) +
  scale_color_brewer(palette="Dark2") +
  coord_flip() +
  labs(y="Consumption (%)",
       x="Height (m)",
       color="Scenario") +
  theme_bw()

ccsp %>%
  filter(strata > 4) %>%
  group_by(run) %>%
  summarize(sum_consump = sum(canopy_consumption_strata))

runs <- c("Fire2-Dry","Fire5-Dry","Fire2-Wet","Fire5-Wet")
spline_list <- list()
for(i in 1:length(runs)){
  spline_int <- as.data.frame(spline(ccsp[ccsp$run==runs[i],]$strata, 
                                     ccsp[ccsp$run==runs[i],]$canopy_consumption_strata))
  spline_int$run <- runs[i]
  spline_list[[i]] <- spline_int
}
spline_df <- bind_rows(spline_list) %>%
  mutate(Rotation = if_else(run %in% c("Fire5-Wet","Fire5-Dry"), "5-Year", "2-Year"),
         Climate = if_else(run %in% c("Fire5-Wet","Fire2-Wet"), "Hot-Wet", "Hot-Dry"))


ccsp_plot <- spline_df %>%
  ggplot() +
  geom_line(aes(x, y, color=Climate, linetype=Rotation)) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("solid","dashed")) +
  scale_x_continuous(limits = c(0,15.3), expand=c(0,0)) +
  scale_y_continuous(limits = c(0,102), expand=c(0,0)) +
  coord_flip() +
  labs(y="Total Fine Fuel Consumption (%)",
       x="Height (m)",
       color="Climate\nScenario",
       linetype="Fire\nRotation") +
  theme_bw()
ccsp_plot
ggsave("consumption_strata.jpg", ccsp_plot, path=here("Plots"),
       height=2.9,width=4.7, units = "in")
