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
(pheno_spp_facet <- pheno %>%
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

(dryas_site_facet <- pheno %>%
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

(sal_site_facet <- pheno %>%
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


(eri_site_facet <- pheno %>%
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

#### CERTAINTY PLOTS ####

(cert_bar_dry <- pheno %>%
   filter(Spp %in% c("DRYINT")) %>%
   filter(obs %in% "phenocam") %>%
   filter(Year %in% 
            c("2019", "2018", "2017", "2016")) %>%
   filter(sex %in% c("F", "M") | is.na(sex)) %>%
   filter(!(cert_ID %in% c("Q1_Cert", 
                           "Q2_Cert", "Q3-M_Cert", "Q3-F_Cert"))) %>%
   filter(!(cert %in% "N") & !is.na(cert)) %>%
   ggplot() +
   aes(x = cert_ID, fill = cert) +
   geom_bar() +
   coord_flip() +
   scale_fill_viridis_d(option = "plasma") +
   labs(x = "Phenophase", y = "Count", fill = "Certainty Index ") +
   theme_classic() +
   facet_wrap(vars(Spp, Year)))

(cert_bar_sal <- pheno %>%
    filter(Spp %in% c("SALARC")) %>%
    filter(obs %in% "phenocam") %>%
    filter(Year %in% 
             c("2019", "2018", "2017", "2016")) %>%
    filter(sex %in% c("F", "M") | is.na(sex)) %>%
    filter(!(cert_ID %in% c("Q1_Cert", 
                            "Q2_Cert", "Q3-M_Cert", "Q3-F_Cert"))) %>%
    filter(!(cert %in% "N") & !is.na(cert)) %>%
    ggplot() +
    aes(x = cert_ID, fill = cert) +
    geom_bar() +
    coord_flip() +
    scale_fill_viridis_d(option = "plasma") +
    labs(x = "Phenophase", y = "Count", fill = "Certainty Index ") +
    theme_classic() +
    facet_wrap(vars(Spp, Year)))
#### ALL DRYAS OBS ####
(dryas_timeline <- pheno %>%
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
(salix_timeline <- pheno %>%
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
(eri_timeline <- pheno %>%
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
  # Filter the data for the correct phase and year range, but keep both 'transect' and 'phenocam'
  filtered_data <- df %>%
    filter(Year >= 2016, Year <= 2019) # First, filter by year

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
               hjust = 1.1, vjust = 1.5, size = 4, color = "black", fontface = "italic")
    
    return(list(plot = plot, p_value = p_value))
  } else {
    # Return an empty plot if there are insufficient levels or no data
    print("Insufficient data or only one level of obs, returning empty plot.")
    return(list(plot = ggplot() + geom_blank(), p_value = NA))
  }
}

# Run the modified function on phases
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


# Save the ANOVA summary table to a CSV file
write.csv(anova_summaries, file = "figures/anova_summary.csv", row.names = FALSE)

#### Overall plot ####
overall <- pheno %>%
  filter(Year >= 2016L & Year <= 2019L) %>%
  filter(phase_ID %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7")) %>%
  filter(Spp %in% c("ERIVAG", "SALARC", "DRYINT","SNOW")) %>%
  drop_na(phase_DATE)

overall$Spp2 <- overall$Spp
overall2 <- overall %>% 
  unite(plot,Spp2, phase_ID)

# boxplot to show distribution 2016-2019
(overall_plot <- overall2 %>%
    ggplot() +
    aes(x = plot, y = phase_DATE, fill = Spp, alpha = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "Phenology stuff", fill = "Observation type") +
    scale_alpha_manual(values=c(1, 0.5)) +
    theme_classic() +
    theme(legend.position = "none"))
overall_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

overall2$plot <- factor (overall2$plot, levels = c("SNOW_P1", "ERIVAG_P1", "ERIVAG_P2", "ERIVAG_P3", "DRYINT_P1", "DRYINT_P2", "DRYINT_P3", "DRYINT_P4", "DRYINT_P5", "DRYINT_P6", "SALARC_P1", "SALARC_P2", "SALARC_P3", "SALARC_P4", "SALARC_P5", "SALARC_P6", "SALARC_P7")) 
(overall_plot <- overall2 %>%
    ggplot() +
    aes(x = phase_DATE, y = plot, fill = Spp, alpha = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    labs(x = "DOY (2016 - 2019)", y = "Phenophase", title = "Phenology stuff", fill = "Observation type") +
    scale_alpha_manual(values=c(1, 0.5)) +
    xlim(0, 365) +
    theme_classic() +
    theme(legend.position = "none"))

(overall_plot_summeronly <- overall2 %>%
    ggplot() +
    aes(x = phase_DATE, y = plot, fill = Spp, alpha = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    labs(x = "DOY (2016 - 2019)", y = "Phenophase", title = "Phenophase Chronology (Transect [light] vs Phenocam [dark])", fill = "Observation type") +
    scale_alpha_manual(values=c(1, 0.5)) +
    xlim(120,290) +
    theme_minimal() +
    theme(legend.position = "none"))

ggsave(overall_plot_summeronly, filename = "figures/QHI_spp_chronto_alt.png", height = 10, width = 12)

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

