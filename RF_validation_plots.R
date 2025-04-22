library(here)
library(tidyverse)
library(ggthemes)
library(tidymodels)

spgrps <- c(1,2,3,4)

df_list <- list()
for(i in 1:length(spgrps)){
  metrics <- read.csv(here("9.FIA","RF_Models","Plots",paste0("cv_metrics_spgrp",spgrps[i],".csv")))
  metrics$spgrp <- factor(spgrps[i])
  df_list[[i]] <- metrics
}

metrics <- bind_rows(df_list)
  
metrics %>% 
  rename(`R-squared` = R2) %>%
  pivot_longer(cols = c("RMSE","R-squared"),
               names_to = "metric",
               values_to = "value") %>%
  mutate(spgrp = factor(spgrp, levels = c(1,2,3,4),
                        labels = c("Pines","Other Conifers","Soft Hardwoods","Hard Hardwoods"))) %>%
  ggplot() +
  geom_boxplot(aes(spgrp,value,color=spgrp)) +
  geom_jitter(aes(spgrp,value,color=spgrp, fill=spgrp), width=0.2, shape=21, alpha=0.5) +
  facet_wrap(.~metric, scales="free_x") +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_x_discrete(limits=rev) +
  labs(x="Species Group",
       y="RMSE") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("model_metrics.png", path=here("9.FIA","RF_Models","Plots"), width=7, height=1.8)

##
df_list <- list()
for(i in 1:length(spgrps)){
  preds <- read.csv(here("9.FIA","RF_Models","Plots",paste0("cv_predictions_spgrp",spgrps[i],".csv")))
  preds$spgrp <- factor(spgrps[i])
  df_list[[i]] <- preds
}
preds <- bind_rows(df_list)

preds <- preds %>%
  mutate(age_class = case_when(y_test < 50 ~ "< 50",
                               y_test < 100 ~ "50-100",
                               y_test < 150 ~ "100-150",
                               y_test < 200 ~ "150-200",
                               TRUE ~ "200+")) %>%
  mutate(age_class = factor(age_class, levels = c("< 50", "50-100","100-150","150-200","200+")),
         spgrp = factor(spgrp, levels = c(1,2,3,4),
                        labels = c("Pines","Other Conifers","Soft Hardwoods","Hard Hardwoods")))

preds %>% 
  ggplot() +
  geom_point(aes(y_test,y_pred,color=factor(Fold)), shape=1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="black") +
  coord_equal() +
  facet_wrap(.~spgrp) +
  labs(y="Predicted Age",
       x="Actual Age",
       color="CV Fold") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("rmse_all.png", path=here("9.FIA","RF_Models","Plots"), height = 5.4, width = 6)

##

rmse_ages <- preds %>%
  group_by(spgrp, age_class, Fold) %>%
  summarize(RMSE = sqrt(mean((y_pred - y_test)^2))) %>%
  ungroup() %>%
  group_by(spgrp, age_class) %>%
  summarize(RMSE_avg = mean(RMSE),
            RMSE_se = sd(RMSE)/sqrt(n()))

rmse_ages %>%
  ggplot() +
  geom_errorbar(aes(age_class, ymin = RMSE_avg-RMSE_se, ymax = RMSE_avg+RMSE_se), width = 0.3) +
  geom_point(aes(age_class, RMSE_avg, color = spgrp, fill = spgrp), shape=21, alpha=0.5) +
  facet_wrap(.~spgrp) +
  labs(x="Age class (years)",
       y="RMSE (years)") +
  scale_y_continuous(limits = c(0,120), 
                     expand=c(0,0),
                     breaks = seq(0,120,20)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  theme_bw() +
  theme(legend.position = "none")

ggsave("rmse_age_classes.png", path=here("9.FIA","RF_Models","Plots"), width = 6, height=4)
