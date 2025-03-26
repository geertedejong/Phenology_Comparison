### 2. Compare transect and phenocam data ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier, Elise Gallois, Anne Bjorkman & Isla Myers-Smith ####
### Date: 26th March 2025 ###

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

#### LOAD FULL PHENOLOGY DATA DATA ####
pheno <- read.csv(file = "data/phenology_transect_cam.csv")

str(pheno)

#### SNOW FREE DAY - S3 to P1 ####
# rename column s3 to p1 to match transect data
pheno <- pheno %>% 
  mutate(pheno, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 

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
   facet_wrap(vars(Spp),nrow = 3,
              ncol = 1,))



#### VISUALISE PHENOPHASES BY SITE ####

(dryas_site_facet <- pheno %>%
   filter(Plot.ID %in% c("1","2","5")) %>%
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
   facet_wrap(vars(Plot.ID)))

(sal_site_facet <- pheno %>%
    filter(Plot.ID %in% c("1","2","3","4","5","6")) %>%
    filter(Spp %in% "SALARC") %>%
    filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6", "P1","P7",
                             'P1_before','P2_before','P3_before','P4_before',
                             'P5_before','P6_before','P7_before'))) %>%
    ggplot() +
    aes(x = phase_ID, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    scale_fill_viridis_d(option = "viridis") +
    theme_bw() +
    labs(x = "Phenophase", title = "Salix Phenophases by Site",y = "DOY (2016-2019)", fill = "Observation type") +
    facet_wrap(vars(Plot.ID)))


(eri_site_facet <- pheno %>%
    filter(Plot.ID %in% c("1","3","6")) %>%
    filter(Spp %in% "ERIVAG") %>%
    filter(!(phase_ID %in% c("S1", "S2", "S3", "S4", "S5", "S6", "P1","P7",
                             'P1_before','P2_before','P3_before','P4_before',
                             'P5_before','P6_before','P7_before'))) %>%
    ggplot() +
    aes(x = phase_ID, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    scale_fill_viridis_d(option = "viridis") +
    theme_bw() +
    labs(x = "Phenophase",title = "Eriophorum Phenophases by Site", y = "DOY (2016-2019)", fill = "Observation type") +
    facet_wrap(vars(Plot.ID)))

#### CERTAINTY PLOTS ####

(cert_bar_dry <- pheno %>%
   filter(Spp %in% c("DRYINT")) %>%
   filter(obs %in% "phenocam") %>%
   
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
   facet_wrap(vars(Spp)))

(cert_bar_sal <- pheno %>%
    filter(Spp %in% c("SALARC")) %>%
    filter(obs %in% "phenocam") %>%
    
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
    facet_wrap(vars(Spp)))
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

#### ANOVA - Frequentist approach ####
# function to run ANOVA and generate boxplots with stats included
anova_boxplot <- function(df, phase_id, species = NULL, title = "") {
  filtered_data <- df %>% 
    filter(phase_ID == phase_id, Year >= 2016 & Year <= 2019) 
  
  if (!is.null(species)) {
    filtered_data <- filtered_data %>% filter(Spp == species)
  }
  
  anova_result <- aov(phase_DATE ~ obs, data = filtered_data)
  anova_summary <- tidy(anova_result)
  
  # extract  f-statistic and p-value
  f_statistic <- anova_summary %>% 
    filter(term == "obs") %>%
    pull(statistic)
  p_value <- anova_summary %>% 
    filter(term == "obs") %>%
    pull(p.value)

  annotation_y <- max(filtered_data$phase_DATE, na.rm = TRUE) + 10
  
  # ggplot function
  plot <- filtered_data %>%
    ggplot(aes(x = obs, y = phase_DATE, fill = obs, col = obs)) +
    geom_boxplot(alpha = 0.8, outlier.colour = NA) + # Add transparency to the boxes
    geom_jitter(width = 0.2, size = 2, alpha = 0.2) + 
    hrbrthemes::scale_fill_ipsum() +
    hrbrthemes::scale_colour_ipsum() +
    #stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", 
         title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none") +
    annotate("text", x = Inf, y = annotation_y, 
             label = paste("F =", round(f_statistic, 2), "\n p =", format.pval(p_value)), 
             hjust = 1.1, vjust = 1.5, size = 4, color = "black", fontface = "italic")
  
  return(list(plot = plot, summary = anova_summary))
}

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

# use function to run analysis for each phenophase and save plots / summaries
results <- map(phases, function(phase) {
  anova_boxplot(pheno, phase$phase_id, phase$species, phase$title)
})

#  plots and summaries
plots <- map(results, "plot")
anova_summaries <- map_dfr(results, "summary")

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

