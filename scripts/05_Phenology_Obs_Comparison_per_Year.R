### 5. Final figures ###
### Geerte Fälthammar de Jong, gugeerte@gmail.com ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ###
### Date: November 2024 ###

#### LOAD PACKAGES  #####

library(dplyr)
library(readr)
library(tidyverse) 
library(esquisse)
library(ggpubr)
library(gridExtra)
library(hrbrthemes)
library(purrr)
library(broom) 
library(lme4)
library(lmerTest)
library(ggpattern)
library(forcats)
library(RColorBrewer)

#### LOAD FULL PHENOLOGY DATA ####
pheno <- read.csv(file = "data/phenology_transect_cam.csv")
pheno_clean <- read.csv(file = "data/phenology_transect_cam_CLEAN.csv")
s2    <- read.csv(file = "data/S2QHIphenocam.csv")

str(pheno_clean)

s2_ndvisf<- subset(s2, NDVI_20m>0.2) #remove all NDVI values below o.2 to exclude negatives and snow
s2_ndsi <- subset(s2, NDSI_20m>0.4)

#### SNOW FREE DAY - S3 to P1 ####
# rename column s3 to p1 to match transect data
pheno <- pheno %>% 
  mutate(pheno, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 


pheno_clean <- pheno_clean %>% 
  mutate(pheno_clean, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 


# filter the dataset to elminate duplicates

# Filtering the dataset
pheno_filtered <- pheno_clean %>%
  filter(!is.na(phase_DATE))             # Keep rows where phase_DATE is not NA

pheno <- pheno %>% select(-Q_ID)
# i want to remove duplicates in pheno
# pheno_clean <- pheno[!duplicated(pheno[c('Spp', 'Plot.ID', 'Year', 'ind.ID', 'phase_ID', 'phase_DATE')]), ]

# list of all phases
phases <- list(
  list(phase_id = "P1", species = NULL, title = "First Day 100% Snow Free"),
  list(phase_id = "P2", species = "ERIVAG", title = "E. vaginatum First Bud Visible"),
  list(phase_id = "P2", species = "DRYINT", title = "D. integrifolia First Bud Visible"),
  list(phase_id = "P3", species = "DRYINT", title = "D. integrifolia First Open Flower"),
  list(phase_id = "P4", species = "DRYINT", title = "D. integrifolia First Petal Shed"),
  list(phase_id = "P5", species = "DRYINT", title = "D. integrifolia First Filament Twist"),
  list(phase_id = "P2", species = "SALARC", title = "S. arctica First Leaf Bud Burst"),
  list(phase_id = "P5", species = "SALARC", title = "S. arctica First Leaf Turns Yellow"),
  list(phase_id = "P6", species = "SALARC", title = "S. arctica Last Leaf Turns Yellow")
)

#### Figure 1 comparative boxplots ####

# Function to generate boxplot and run model
anova_boxplot <- function(df, phase_id, species = NULL, title = "") {
  # Filter data
  filtered_data <- df %>%
    filter(phase_ID == phase_id, Year %in% 2016:2019, if (!is.null(species)) Spp == species else TRUE) %>%
    mutate(obs = ifelse(obs == "transect", "in-situ", obs))
  
  # Calculate summary statistics
  doy_stats <- filtered_data %>%
    group_by(obs) %>%
    summarise(mean_DOY = mean(phase_DATE, na.rm = TRUE),
              sd_DOY = sd(phase_DATE, na.rm = TRUE), .groups = "drop")
  
  mean_diff <- abs(diff(doy_stats$mean_DOY))
  sd_diff <- sqrt(mean(doy_stats$sd_DOY^2))
  
  # Fit mixed models
  full_model <- lmer(phase_DATE ~ obs + (1|Year), data = filtered_data, na.action = na.omit)
  null_model <- lmer(phase_DATE ~ 1 + (1|Year), data = filtered_data, na.action = na.omit)
  
  # Likelihood ratio test
  lrt <- anova(full_model, null_model)
  p_value <- lrt[2, "Pr(>Chisq)"]
  
  # Model predictions
  newdat <- expand.grid(obs = unique(filtered_data$obs)) %>%
    mutate(phase_DATE = predict(full_model, ., re.form = NA))
  
  mm <- model.matrix(terms(full_model), newdat)
  pvarl <- diag(mm %*% tcrossprod(vcov(full_model), mm))
  
  newdat <- newdat %>%
    mutate(plo = phase_DATE - 2 * sqrt(pvarl),
           phi = phase_DATE + 2 * sqrt(pvarl))
  
  # Create plot
  plot <- ggplot(filtered_data, aes(x = obs, y = phase_DATE, fill = obs, col = obs)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.2) +
    geom_pointrange(data = newdat, aes(x = obs, y = phase_DATE, ymin = plo, ymax = phi, color = obs, fill = obs),
                    size = 1, fatten = 1.5, inherit.aes = FALSE) +
    hrbrthemes::scale_fill_ipsum() + hrbrthemes::scale_colour_ipsum() +
    labs(x = "Observation method", y = "DOY (2016 - 2019)", title = title) +
    theme_classic() + theme(legend.position = "right") +
    annotate("text", x = Inf, y = max(filtered_data$phase_DATE, na.rm = TRUE),
             label = paste("p =", format.pval(p_value)), hjust = 1.1, vjust = 0.3, size = 4, fontface = "italic")
  
  return(list(plot = plot, p_value = p_value, mean_diff = mean_diff, sd_diff = sd_diff))
}

anova_boxplot <- function(df, phase_id, species = NULL, title = "") {
  # Filter data
  filtered_data <- df %>%
    filter(phase_ID == phase_id, Year %in% 2016:2019, if (!is.null(species)) Spp == species else TRUE) %>%
    mutate(obs = ifelse(obs == "transect", "in-situ", obs))
  
  # Calculate summary statistics
  doy_stats <- filtered_data %>%
    group_by(obs) %>%
    summarise(mean_DOY = mean(phase_DATE, na.rm = TRUE),
              sd_DOY = sd(phase_DATE, na.rm = TRUE), .groups = "drop")
  
  mean_diff <- abs(diff(doy_stats$mean_DOY))
  sd_diff <- sqrt(mean(doy_stats$sd_DOY^2))
  
  # Fit mixed models
  full_model <- lmer(phase_DATE ~ obs + (1|Year), data = filtered_data, na.action = na.omit)
  null_model <- lmer(phase_DATE ~ 1 + (1|Year), data = filtered_data, na.action = na.omit)
  
  # Likelihood ratio test
  lrt <- anova(full_model, null_model)
  p_value <- lrt[2, "Pr(>Chisq)"]
  
  # Model predictions
  newdat <- expand.grid(obs = unique(filtered_data$obs)) %>%
    mutate(phase_DATE = predict(full_model, ., re.form = NA))
  
  mm <- model.matrix(terms(full_model), newdat)
  pvarl <- diag(mm %*% tcrossprod(vcov(full_model), mm))
  
  newdat <- newdat %>%
    mutate(plo = phase_DATE - 2 * sqrt(pvarl),
           phi = phase_DATE + 2 * sqrt(pvarl))
  
  # Create plot
  plot <- ggplot(filtered_data, aes(x = obs, y = phase_DATE, fill = obs, col = obs)) +
    geom_jitter(width = 0.2, size = 4, alpha = 0.3, aes(shape = "Observed")) +
    geom_pointrange(data = newdat, aes(x = obs, y = phase_DATE, ymin = plo, ymax = phi, group = obs, fill = obs, shape = "Modeled"), colour = "black",
                    size = 4, fatten = 2, inherit.aes = FALSE) +
    geom_pointrange(data = newdat, aes(x = obs, y = phase_DATE, ymin = plo, ymax = phi, color = obs, fill = obs, shape = "Modeled"),
                    size = 3, fatten = 2, inherit.aes = FALSE) +
    scale_shape_manual(name = "Data Type", values = c("Observed" = 16, "Modeled" = 17)) +
    hrbrthemes::scale_fill_ipsum() + hrbrthemes::scale_colour_ipsum() +
    labs(x = "", y = "DOY (2016 - 2019)", title = title, shape = "Data Type") +
    theme_classic() + theme(legend.position = "none") +
    annotate("text", x = Inf, y = max(filtered_data$phase_DATE, na.rm = TRUE), label = paste("p =", format.pval(round(p_value, digits=3),digits=3, nsmall = 3)), hjust = 1.1, vjust = 0.3, size = 6, fontface = "italic")+
    theme(plot.title=element_text(color="black", size=18, vjust=1.25)) +
    theme(axis.text.x=element_text(size=20,color="black", margin=margin(4,4,5,4,"pt"))) +
    theme(axis.text.y=element_text(size=20,color="black", margin=margin(4,4,5,4,"pt"))) +
    theme(axis.title.x=element_text(size=20,color="black", margin=margin(10,0,0,0))) +
    theme(axis.title.y=element_text(size=22,color="black", margin=margin(0,10,0,0)))
  
  return(list(plot = plot, p_value = p_value, mean_diff = mean_diff, sd_diff = sd_diff))
}



# Run the function on phases
results <- map(phases, function(phase) {
  anova_boxplot(pheno_clean, phase$phase_id, phase$species, phase$title)
})

# Extract plots, p-values, and summaries
plots <- map(results, "plot")
p_values <- map(results, "p_value")
mean_diffs <- map(results, "mean_diff")
sd_diffs <- map(results, "sd_diff")

# Arrange the plots into a grid
anova_pheno <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)

anova_pheno

anova_pheno_long <- ggarrange(plotlist = plots, ncol = 2, nrow = 5)

# save the full panel of ANOVA plots
ggsave(anova_pheno, filename = "figures/anova_phenocam_box_2024.png", height = 13, width = 15)

ggsave(anova_pheno_long, filename = "figures/anova_phenocam_box_2024_long.png", height = 18, width = 10)


#### Figure 3 Boxplots chronologically ####
overall <- pheno_clean %>%
  filter(Year >= 2016L & Year <= 2019L) %>%
  filter(phase_ID %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7")) %>%
  filter(Spp %in% c("ERIVAG", "SALARC", "DRYINT","SNOW")) %>%
  drop_na(phase_DATE)

overall$Spp2 <- overall$Spp
overall2 <- overall %>% 
  unite(plot,Spp2, phase_ID)

# boxplot to show distribution 2016-2019
# Set factor levels
overall2$plot <- factor(overall2$plot, 
                        levels = c("SNOW_P1", "ERIVAG_P1", "ERIVAG_P2", "ERIVAG_P3", "DRYINT_P1", "DRYINT_P2", 
                                   "DRYINT_P3", "DRYINT_P4", "DRYINT_P5", "DRYINT_P6", "SALARC_P1", "SALARC_P2", 
                                   "SALARC_P3", "SALARC_P4", "SALARC_P5", "SALARC_P6", "SALARC_P7"))

# Rename levels using fct_recode
overall2$plot <- fct_recode(overall2$plot,
                            "Snow Free Plot" = "SNOW_P1",
                            "E.vaginatum Snow Free" = "ERIVAG_P1",
                            "E.vaginatum First Flower Bud" = "ERIVAG_P2",
                            "E.vaginatum First Pollen visible" = "ERIVAG_P3",
                            "D.integrifolia Snow Free" = "DRYINT_P1",
                            "D.integrifolia First White on Bud" = "DRYINT_P2",
                            "D.integrifolia First Open Flower" = "DRYINT_P3",
                            "D.integrifolia Last Petal Shed" = "DRYINT_P4",
                            "D.integrifolia First Twisting of Filaments" = "DRYINT_P5",
                            "D.integrifolia All leaves dead" = "DRYINT_P6",
                            "S.arctica Snow Free" = "SALARC_P1",
                            "S.arctica First Leaf Bud Burst " = "SALARC_P2",
                            "S.arctica First Pollen" = "SALARC_P3",
                            "S.arctica Onset Seed Dispersal" = "SALARC_P4",
                            "S.arctica First Yellow Leaf" = "SALARC_P5",
                            "S.arctica Last Green Leaf" = "SALARC_P6",
                            "S.arctica All Leaves Dead" = "SALARC_P7")


plot_boxplot <- function(data, x_var, y_var, x_limits = NULL, title = "Phenology stuff", xlabel = "DOY (2016 - 2019)", ylabel = "Phenophase") {
  data %>%
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }}, fill = interaction(Spp, obs))) +
    geom_boxplot() +
    #scale_fill_brewer(palette = "Set2") + # Use a qualitative palette that adapts to number of groups
    scale_fill_manual(values = c("lightgreen","darkgreen","yellow","orange","pink","purple","blue"),
                      breaks = c("SALARC.transect","SALARC.phenocam","DRYINT.transect","DRYINT.phenocam"  ,"ERIVAG.transect","ERIVAG.phenocam","SNOW.phenocam"),
                      labels = c("S.arctica transect","S.arctica phenocam","D.integrifolia transect","D.integrifolia  phenocam","E.vaginatum transect","E.vaginatum phenocam","Snow phenocam"))+
    labs(x = xlabel, y = ylabel, title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "right") +
    { if (!is.null(x_limits)) xlim(x_limits) else NULL }
}

# Full date range plot
overall_plot <- plot_boxplot(overall2, phase_DATE, plot)
overall_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Summer-only plot
overall_plot_summeronly <- plot_boxplot(overall2, phase_DATE, plot, x_limits = c(120, 290), 
                                        title = "Phenophase Chronology")
overall_plot_summeronly

#### Figure 2 Satellites and snow-free obs ####

s2    <- read.csv(file = "data/S2QHIphenocam.csv")

## prep dash lines
# filter for snow free on camera
snowfree16 <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year== 2016L) 
snowfree17 <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year== 2017L)
snowfree18 <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year== 2018L)
snowfree19 <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year== 2019L)


# filter for senescence
salsen2_16 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P6") %>%
  filter(Year ==2016L) 
salsen2_17 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P6") %>%
  filter(Year ==2017L) 
salsen2_18 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P6") %>%
  filter(Year ==2018L) 
salsen1_19 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P5") %>%
  filter(Year ==2019L) 


#### some cleaning ####
s2_ndvisf<- subset(s2, NDVI_20m>0.2) #remove all NDVI values below o.1 to exclude negatives and snow

#### Exploration plot of S2 data ####
# NDVI
# average geom smooth and colors for camera locations, separate years
(s2_plot <- s2_ndvisf %>%
   ggplot() +
   aes(x = doi, y = NDVI_20m) +
   geom_smooth() +
   geom_point(aes(color=name))+
   hrbrthemes::scale_fill_ipsum() +
   labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
   xlim(120,290) +
   facet_grid(year~.) +
   theme_classic() +
   theme(legend.position ="right"))

# average geom smooth and colors for camera locations, average for all years
(s2_plot <- s2_ndvisf %>%
    ggplot() +
    aes(x = doi, y = NDVI_20m) +
    geom_smooth() +
    geom_point(aes(color=year))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
    xlim(120,290) +
    scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "NDVI")) +
    theme_classic() +
    theme(legend.position = "none"))

# NDSI
# average geom smooth and colors for camera locations, separate years
(s2_plot <- s2_ndsi %>%
    ggplot() +
    aes(x = doi, y = NDSI_20m) +
    geom_smooth() +
    geom_point(aes(color=name))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDSI", x = "DOY (2016 - 2019)", fill = "year") +
    facet_grid(year~.) +
    theme_classic() +
    theme(legend.position = "none"))

# average geom smooth and colors for camera locations, average for all years
(s2_plot <- s2_ndsi %>%
    ggplot() +
    aes(x = doi, y = NDSI_20m) +
    geom_smooth() +
    geom_point(aes(color=year))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDSI", x = "DOY (2016 - 2019)", fill = "year") +
    theme_classic() +
    theme(legend.position = "none"))

# NB run 02_Phenology_Obs_Comparison script first to load data into environment
(cam_sf16 <- min(snowfree16$phase_DATE, na.rm=TRUE))
(cam_sf17 <- min(snowfree17$phase_DATE, na.rm=TRUE))
(cam_sf18 <- min(snowfree18$phase_DATE, na.rm=TRUE))
(cam_sf19 <- min(snowfree19$phase_DATE, na.rm=TRUE))
(cam_senescence16 <- max(salsen2_16$phase_DATE, na.rm=TRUE))
(cam_senescence17 <- max(salsen2_17$phase_DATE, na.rm=TRUE))
(cam_senescence18 <- max(salsen2_18$phase_DATE, na.rm=TRUE))
(cam_senescence19 <- max(salsen1_19$phase_DATE, na.rm=TRUE))

# for first doy of ndvi>0.2 per year I manually checked the table

#### combination plot of cams, obs and NDVI for S2####
(comb_plot <- ggplot()+
   #geom_point(data=s2_ndvisf, aes(x=doi, y=NDVI_20m, color=factor(year)), alpha=0.3, size=1)+
   geom_smooth(data=s2_ndvisf,aes(x=doi, y=NDVI_20m, color=factor(year), fill=factor(year)), alpha=0.2) +
   hrbrthemes::scale_color_ipsum() +geom_vline(xintercept=cam_sf16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_senescence16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_sf17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_senescence17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_sf18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_senescence18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_sf19, linetype='dashed',color='purple',size=1)+
   xlim(100,300)+
   ylim(0.1,1)+
   labs(y = "NDVI", x = "DOY (2016 - 2019)", color= "year") +
   annotate(x = cam_sf16-20, y = 1, label = "snow-free cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_sf16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf19, y = 0.76, label = "2019", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence16-20, y = 1, label = "senescence cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_senescence16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   theme_classic() +
   scale_color_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   scale_fill_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   theme(legend.position = "right")
)

ggsave(comb_plot, filename = "figures/cam_obs_ndvi_greencurv.png", height = 10, width = 12)

#### combination plot of cams, obs and NDSI for S2####
(comb_plot <- ggplot()+
   #geom_point(data=s2_ndvisf, aes(x=doi, y=NDVI_20m, color=factor(year)), alpha=0.3, size=1)+
   geom_smooth(data=s2_ndsi,aes(x=doi, y=NDSI_20m, color=factor(year),fill=factor(year)), alpha=0.2) +
   hrbrthemes::scale_color_ipsum() +
   geom_vline(xintercept=cam_sf16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_senescence16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_sf17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_senescence17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_sf18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_senescence18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_sf19, linetype='dashed',color='purple',size=1)+
   xlim(100,300)+
   ylim(0.1,1)+
   labs(y = "NDSI", x = "DOY (2016 - 2019)", color= "year") +
   annotate(x = cam_sf16-20, y = 1, label = "snow-free cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_sf16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf19, y = 0.76, label = "2019", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence16-20, y = 1, label = "senescence cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_senescence16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   theme_classic() +
   scale_color_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   scale_fill_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   theme(legend.position = "right")
)

ggsave(comb_plot, filename = "figures/cam_obs_ndsi_greencurv.png", height = 10, width = 12)

#### combination plot of cams, obs and NDVI for modis####
(comb_plot <- ggplot()+
   #geom_point(data=s2_ndvisf, aes(x=doi, y=NDVI_20m, color=factor(year)), alpha=0.3, size=1)+
   geom_smooth(data=ndvi_m,aes(x=doy, y=NDVI, color=factor(year),fill=factor(year)), alpha=0.2) +
   hrbrthemes::scale_color_ipsum() +
   geom_vline(xintercept=cam_sf16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_senescence16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_sf17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_senescence17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_sf18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_senescence18, linetype='dashed',color='blue',size=1)+
   geom_vline(xintercept=cam_sf19, linetype='dashed',color='purple',size=1)+
   xlim(100,300)+
   ylim(0.1,1)+
   labs(y = "NDVI", x = "DOY (2016 - 2019)", color= "year") +
   annotate(x = cam_sf16-20, y = 1, label = "snow-free cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_sf16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_sf19, y = 0.76, label = "2019", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence16-20, y = 1, label = "senescence cam", vjust = 0, geom = "label",size = 3)+
   annotate(x = cam_senescence16, y = 1, label = "2016", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence17, y = 0.92, label = "2017", vjust = 0, angle= 90, geom = "label",size = 3)+
   annotate(x = cam_senescence18, y = 0.84, label = "2018", vjust = 0, angle= 90, geom = "label",size = 3)+
   theme_classic() +
   scale_color_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   scale_fill_manual(name= "Year", labels = c("2016", "2017", "2018", "2019"), values= c("orange", "green", "blue", "purple"))+
   theme(legend.position = "right")
)

ggsave(comb_plot, filename = "figures/cam_obs_ndvi_greencurv_modis.png", height = 10, width = 12)

#### Figure 3 big plot combining phenocams, transect and sat data ####

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(patchwork)

phenophase_boxplot <- overall2 %>%
    ggplot(aes(x = phase_DATE, y = plot)) +
    geom_boxplot(aes(fill = interaction(Spp, obs))) +
    #geom_smooth(data = s2_ndvisf,aes(x = doi, y = (NDVI_20m)*14, color = factor(year))) +
    #scale_y_continuous(sec.axis = ~(.*14), name = "NDVI_20m")+
    scale_fill_manual(values = c("lightgreen","darkgreen","yellow","orange","pink","purple","blue"),
                      breaks = c("SALARC.transect","SALARC.phenocam","DRYINT.transect","DRYINT.phenocam",
                                 "ERIVAG.transect","ERIVAG.phenocam","SNOW.phenocam"),
                      labels = c("Sal.Arc. transect","Sal.Arc. phenocam","Dry.Int. transect",
                                 "Dry.Int. phenocam","Eri.Vag. transect","Eri.Vag phenocam",
                                 "Snow phenocam")) +
    labs(x = "DOY",  fill = "Observation type") +
    theme_classic()+
    coord_cartesian(xlim = c(130, 280))+
    theme(legend.position = "right")
  


# Function to create the NDVI line plot
plot_ndvi_line <- function(s2_ndvisf, x_limits = NULL, xlabel = "DOY", ylabel = "NDVI") {
  p <- s2_ndvisf %>%
    ggplot() +
    geom_smooth(aes(x = doi, y = NDVI_20m, color = factor(year))) +
    scale_color_ipsum() +
    labs(x = xlabel, y = ylabel, color = "Year") +
    theme_classic() +
    theme(legend.position = "right") +
    coord_cartesian(xlim = x_limits)  # Ensure x-axis matches phenophase plot
  
  return(p)
}

# Define the x-axis limits (e.g., full year or summer only)
x_limits <- c(130, 280)  # For summer-only

# Create both plots with matching x-axis limits
#phenophase_plot <- plot_phenophase_boxplot(overall2, phase_DATE, plot, x_limits = x_limits, title = "Phenophase Chronology")
ndvi_plot <- plot_ndvi_line(s2_ndvisf, x_limits = x_limits)

# Combine plots vertically using patchwork
combined_plot <- phenophase_boxplot / ndvi_plot

# Display the combined plot
combined_plot


#### test to combine plots in one with a secondary axis ####

# Create the primary NDVI line plot
ndvi_plot <- ggplot(s2_ndvisf, aes(x = doi)) +
  geom_smooth(aes(y = NDVI_20m, color = factor(year))) +
  scale_y_continuous(name = "NDVI",
                     sec.axis = sec_axis(~ ., name = "Phenophase")) + # Secondary axis
  scale_color_manual(values = c("blue", "green", "orange", "purple")) +
  theme_classic() +
  labs(color = "Year")

# Create the secondary phenophase boxplot and flip it horizontally
phenophase_plot <- ggplot(overall2, aes(x = phase_DATE, y = plot)) +
  geom_boxplot(aes(fill = interaction(Spp, obs))) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "yellow", "orange", "pink", "purple", "blue"),
                    breaks = c("SALARC.transect", "SALARC.phenocam", "DRYINT.transect", "DRYINT.phenocam",
                               "ERIVAG.transect", "ERIVAG.phenocam", "SNOW.phenocam"),
                    labels = c("Sal.Arc. transect", "Sal.Arc. phenocam", "Dry.Int. transect",
                               "Dry.Int. phenocam", "Eri.Vag. transect", "Eri.Vag phenocam",
                               "Snow phenocam")) +
  theme_classic() +
  labs(fill = "Observation type")

# Combine the plots using ggplot's layering
combined_plot <-phenophase_plot + ndvi_plot
  

# Display the combined plot
print(combined_plot)
# prints them side by side
# new test

# Create the primary NDVI line plot


# Create the secondary phenophase boxplot and flip it horizontally
phenophase_plot <- ggplot(overall2, aes(x = phase_DATE, y = plot)) +
  geom_boxplot(aes(fill = interaction(Spp, obs))) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "yellow", "orange", "pink", "purple", "blue"),
                    breaks = c("SALARC.transect", "SALARC.phenocam", "DRYINT.transect", "DRYINT.phenocam",
                               "ERIVAG.transect", "ERIVAG.phenocam", "SNOW.phenocam"),
                    labels = c("Sal.Arc. transect", "Sal.Arc. phenocam", "Dry.Int. transect",
                               "Dry.Int. phenocam", "Eri.Vag. transect", "Eri.Vag phenocam",
                               "Snow phenocam")) +
  theme_classic() +
  labs(fill = "Observation type")+
  geom_smooth(data=s2_ndvisf,aes(x=doi, y = NDVI_20m, color = factor(year))) +
  scale_y_continuous(name = "NDVI",
                     sec.axis = sec_axis(~ ., name = "Phenophase")) + # Secondary axis
  scale_color_manual(values = c("blue", "green", "orange", "purple")) +
  theme_classic() +
  labs(color = "Year")

print(phenophase_plot)

# Combine the plots using ggplot's layering
combined_plot <-phenophase_plot + ndvi_plot


# Display the combined plot
print(combined_plot)

#### calculate recurrence times for sentinel-2 and modis
# sentinel
#transform doi and year columns into dates
s2_ndvisf$date <- as.Date(s2_ndvisf$doi-1, origin=paste0(s2_ndvisf$year,"-01-01"))
#calculate difference in days between consecutive days
s2_ndvisf$recurrence <- c(NA,diff(s2_ndvisf$date))
#calculate min, max and average, ignoring 0's as those are from different tiles taken the same day
s2_recurrence <- na.omit(s2_ndvisf$recurrence)
s2_recurrence <- s2_recurrence[s2_recurrence!=0]

min_rec <- min(s2_recurrence)
max_rec <- max(s2_recurrence)
avg_rec <- mean(s2_recurrence)
sd_rec <- sd(s2_recurrence)

s2_rec_summer <- s2_recurrence[s2_recurrence <200]
min_sum <- min(s2_rec_summer)
max_sum <- max(s2_rec_summer)
avg_sum <- mean(s2_rec_summer)
sd_sum <- sd(s2_rec_summer)
