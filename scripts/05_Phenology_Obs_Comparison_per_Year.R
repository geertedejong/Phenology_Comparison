### 5. Visualize Phenology_Obs comparison per year ###
### Geerte Fälthammar de Jong, gugeerte@gmail.com ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ###
### Date: September 2024 ###

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

#### LOAD FULL PHENOLOGY DATA DATA ####
pheno <- read.csv(file = "data/phenology_transect_cam.csv")
pheno_clean <- read.csv(file = "data/phenology_transect_cam_CLEAN.csv")

str(pheno_clean)

#### SNOW FREE DAY - S3 to P1 ####
# rename column s3 to p1 to match transect data
pheno <- pheno %>% 
  mutate(pheno, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 


pheno_clean <- pheno_clean %>% 
  mutate(pheno_clean, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 


#### VISUALISE PHENOPHASES BY OBSERVATION TYPE ####
#
(pheno_spp_facet <- pheno_clean %>%
   filter(Spp %in% c("DRYINT","SALARC","ERIVAG")) %>%
   filter(Year %in% 
            c("2019", "2018", "2017", "2016")) %>%
   filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6",
                            'P1_before','P2_before','P3_before','P4_before',
                            'P5_before','P6_before','P7_before'))) %>%
   ggplot() +
   aes(x = phase_ID, y = phase_DATE, fill = obs) +
   geom_boxplot() +
   scale_fill_viridis_d(option = "viridis") +
   labs(x = "Phenophase", y = "DOY (2016-2019)", fill = "Observation type") +
   theme_bw() +
   facet_wrap(vars(Spp, Year),nrow = 4,
              ncol = 3,))



#### VISUALISE PHENOPHASES BY SITE ####

(dryas_site_facet <- pheno_clean %>%
   filter(Plot.ID %in% c("1","2","5")) %>%
   filter(Year %in% 
            c("2019", "2018", "2017", "2016")) %>%
   filter(Spp %in% "DRYINT") %>%
   filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6", "P1","P7",
                            'P1_before','P2_before','P3_before','P4_before',
                            'P5_before','P6_before','P7_before'))) %>%
   ggplot() +
   aes(x = phase_ID, y = phase_DATE, fill = obs) +
   geom_boxplot() +
   scale_fill_viridis_d(option = "viridis") +
   theme_bw() +
   labs(x = "Phenophase",title = "Dryas Phenophases by Site", y = "DOY (2016-2019)", fill = "Observation type") +
   facet_wrap(vars(Plot.ID, Year)))

(sal_site_facet <- pheno_clean %>%
    filter(Plot.ID %in% c("1","2","3","4","5","6")) %>%
    filter(Spp %in% "SALARC") %>%
    filter(Year %in% 
             c("2019", "2018", "2017", "2016")) %>%
    filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6", "P1","P7",
                             'P1_before','P2_before','P3_before','P4_before',
                             'P5_before','P6_before','P7_before'))) %>%
    ggplot() +
    aes(x = phase_ID, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    scale_fill_viridis_d(option = "viridis") +
    theme_bw() +
    labs(x = "Phenophase", title = "Salix Phenophases by Site",y = "DOY (2016-2019)", fill = "Observation type") +
    facet_wrap(vars(Plot.ID, Year)))


(eri_site_facet <- pheno_clean %>%
    filter(Plot.ID %in% c("1","3","6")) %>%
    filter(Spp %in% "ERIVAG") %>%
    filter(Year %in% 
             c("2019", "2018", "2017", "2016")) %>%
    filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6", "P1","P7",
                             'P1_before','P2_before','P3_before','P4_before',
                             'P5_before','P6_before','P7_before'))) %>%
    ggplot() +
    aes(x = phase_ID, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    scale_fill_viridis_d(option = "viridis") +
    theme_bw() +
    labs(x = "Phenophase",title = "Eriophorum Phenophases by Site", y = "DOY (2016-2019)", fill = "Observation type") +
    facet_wrap(vars(Plot.ID, Year)))


#### ALL DRYAS OBS ####
(dryas_timeline <- pheno_clean %>%
   filter(Spp %in% "DRYINT") %>%
   filter(phase_ID %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7")) %>%
   ggplot() +
   aes(x = Year, y = phase_DATE, colour = phase_ID) +
   geom_point(size = 1L) +
   stat_smooth(method = 'lm') +
   scale_color_brewer(palette = "Set2") +
   labs(y = "DOY", title = "Dryas Phenophase Shifts (2016 - 2019)", col = "Phenophase") +
   theme_bw())


#### ALL SALIX OBS ####
(salix_timeline <- pheno_clean %>%
   filter(Spp %in% "SALARC") %>%
   filter(phase_ID %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7")) %>%
   ggplot() +
   aes(x = Year, y = phase_DATE, colour = phase_ID) +
   geom_point(size = 1L) +
   stat_smooth(method = 'lm') +
   scale_color_brewer(palette = "Set2") +
   labs(y = "DOY", title = "Salix Phenophase Shifts (2001 - 2019)", col = "Phenophase") +
   theme_bw())

#### ALL ERIOPHORUM OBS ####
(eri_timeline <- pheno_clean %>%
   filter(Spp %in% "ERIVAG") %>%
   filter(phase_ID %in% c("P1", "P2", "P3")) %>%
   ggplot() +
   aes(x = Year, y = phase_DATE, colour = phase_ID) +
   geom_point(size = 1L) +
   stat_smooth(method = 'lm') +
   scale_color_brewer(palette = "Set2") +
   labs(y = "DOY", title = "Eriophorum Phenophase Shifts (2001 - 2019)", col = "Phenophase") +
   theme_bw())

# filter the dataset to elminate duplicates

# Filtering the dataset
pheno_filtered <- pheno_clean %>%
  filter(!is.na(phase_DATE))             # Keep rows where phase_DATE is not NA

pheno <- pheno %>% select(-Q_ID)
# i want to remove duplicates in pheno
pheno_clean <- pheno[!duplicated(pheno[c('Spp', 'Plot.ID', 'Year', 'ind.ID', 'phase_ID', 'phase_DATE')]), ]

# list of all phases
phases <- list(
  list(phase_id = "P1", species = NULL, title = "First Day 100% Snow Free"),
  list(phase_id = "P2", species = "ERIVAG", title = "First E. vaginatum Bud Appearance"),
  list(phase_id = "P2", species = "DRYINT", title = "First D. integrifolia Bud Appearance"),
  list(phase_id = "P3", species = "DRYINT", title = "First D. integrifolia Open Flower"),
  list(phase_id = "P4", species = "DRYINT", title = "First D. integrifolia Petal Shed"),
  list(phase_id = "P5", species = "DRYINT", title = "First D. integrifolia Twisting of Filament"),
  list(phase_id = "P2", species = "SALARC", title = "S. arctica First Leaf Bud Burst"),
  list(phase_id = "P5", species = "SALARC", title = "S. arctica First Leaf Turns Yellow"),
  list(phase_id = "P6", species = "SALARC", title = "S. arctica Last Leaf Turns Yellow")
)

# Function to generate boxplot and run model
anova_boxplot <- function(df, phase_id, species = NULL, title = "") {
  # Filter the data for the correct phase and year range, and filter by species
  filtered_data <- df %>%
    filter(phase_ID == phase_id,  # Filter by phase ID
           Year >= 2016, Year <= 2019,  # Filter by year range
           if (!is.null(species)) Spp == species else TRUE)  # filter by species
  
  # Debugging: Print filtered data summary
  print(paste("Phase:", phase_id, "Species:", ifelse(is.null(species), "All", species)))
  print(paste("Number of rows:", nrow(filtered_data)))
  print(paste("Unique obs:", unique(filtered_data$obs)))
  
  # Check if the filtered data is empty or obs has less than 2 levels
  if (nrow(filtered_data) > 0 && length(unique(filtered_data$obs)) > 1) {
    # Use lmer for mixed model
    full_model <- lmer(phase_DATE ~ obs + (1|Year), data = filtered_data, na.action = na.omit)
    null_model <- lmer(phase_DATE ~ 1 + (1|Year), data = filtered_data, na.action = na.omit)
    
    # Likelihood ratio test
    lrt_result <- anova(full_model, null_model)
    chi_square_statistic <- lrt_result[2, "Chisq"]
    p_value <- lrt_result[2, "Pr(>Chisq)"]
    
    annotation_y <- max(filtered_data$phase_DATE, na.rm = TRUE)
    
    # ggplot function
    plot <- filtered_data %>%
      ggplot(aes(x = obs, y = phase_DATE, fill = obs, col = obs)) +
      geom_boxplot(alpha = 0.8, outlier.colour = NA) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.2) + 
      hrbrthemes::scale_fill_ipsum() +
      hrbrthemes::scale_colour_ipsum() +
      labs(x = "Observation type", y = "DOY (2016 - 2019)", 
           title = title, fill = "Observation type") +
      theme_classic() +
      theme(legend.position = "none") +
      annotate("text", x = Inf, y = annotation_y, 
               label = paste("\n p =", format.pval(p_value)), 
               hjust = 1.1, vjust = 0.3, size = 4, color = "black", fontface = "italic")
    
    return(list(plot = plot, p_value = p_value))
  } else {
    # Return an empty plot if there are insufficient levels or no data
    print("Insufficient data or only one level of obs, returning empty plot.")
    return(list(plot = ggplot() + geom_blank(), p_value = NA))
  }
}


# Run the function on phases
results <- map(phases, function(phase) {
  anova_boxplot(pheno_clean, phase$phase_id, phase$species, phase$title)
})


# Extract plots and summaries
plots <- map(results, "plot")
p_value <- map(results, "p_value")

# Arrange the plots into a grid
anova_pheno <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)

anova_pheno



# save the full panel of ANOVA plots
ggsave(anova_pheno, filename = "figures/anova_phenocam_box_2024.png", height = 10, width = 12)


#### Overall plot ####
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

# General plot function
plot_boxplot <- function(data, x_var, y_var, x_limits = NULL, title = "Phenology stuff", xlabel = "DOY (2016 - 2019)", ylabel = "Phenophase") {
  data %>%
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }}, fill = Spp, alpha = obs)) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    scale_alpha_manual(values = c(1, 0.5)) +
    labs(x = xlabel, y = ylabel, title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none") +
    { if (!is.null(x_limits)) xlim(x_limits) else NULL }
}

# Full date range plot
overall_plot <- plot_boxplot(overall2, phase_DATE, plot)
overall_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Summer-only plot
overall_plot_summeronly <- plot_boxplot(overall2, phase_DATE, plot, x_limits = c(120, 290), 
                                        title = "Phenophase Chronology (Transect [light] vs Phenocam [dark])")

# Save summer-only plot
ggsave(overall_plot_summeronly, filename = "figures/QHI_spp_chronto_alt.png", height = 10, width = 12)


### test
library(RColorBrewer) # Load the library for color palettes

plot_boxplot <- function(data, x_var, y_var, x_limits = NULL, title = "Phenology stuff", xlabel = "DOY (2016 - 2019)", ylabel = "Phenophase") {
  data %>%
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }}, fill = interaction(Spp, obs))) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set2") + # Use a qualitative palette that adapts to number of groups
    scale_alpha_manual(values = c(1, 0.5)) +
    labs(x = xlabel, y = ylabel, title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "right") +
    { if (!is.null(x_limits)) xlim(x_limits) else NULL }
}


### end test



#### all pheno phases chrono ###

(fulltimeline_chrono <- pheno %>%
    filter(!(Spp %in% c("SNOW",  "LUPARC", "PEDIC"))) %>%
    filter(!(phase_ID %in% c("S1", "S4","S5","S6","S2", "P1_before", "P2_before", "P3_before", "P4_before", 
                             "P5_before", "P6_before", "P7_before"))) %>%
    ggplot() +
    aes(x = Year, y = phase_DATE, colour = phase_ID) +
    geom_point(size = 1L) +
    stat_smooth(method = 'lm') +
    theme_minimal() +
    labs(x = "Year", y = "DOY", fill = "Phenophase ID") +
    
    facet_wrap(vars(Spp)))

ggsave(fulltimeline_chrono, filename = "figures/QHI_spp_chronto.png", height = 10, width = 12)

