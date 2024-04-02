### 2. Compare transect and phenocam data ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ####
### Date: 17th February 2021 ###

#### LOAD PACKAGES  #####

library(dplyr)
library(readr)
library(tidyverse) 
library(esquisse)
library(ggpubr)
library(gridExtra)

#### LOAD FULL PHENOLOGY DATA DATA ####
pheno <- read.csv(file = "data/phenology/phenology_transect_cam.csv")

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

# P1 - Snow Free Date

# filter for only P1
snowfree <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
snowfree_anova <- aov(phase_DATE ~ obs, data = snowfree)
summary(snowfree_anova)

# boxplot to show distribution 2016-2019
(snow_free <- snowfree %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First Day 100% Snow Free", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Eriophorum P2 - Appearance of first flower bud

# filter for only P2
eriobud <- pheno %>% filter(Spp %in% "ERIVAG") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
eriobud_anova <- aov(phase_DATE ~ obs, data = eriobud)
summary(eriobud_anova)

# boxplot to show distribution 2016-2019
(eriobud_plot <- eriobud %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First E. vaginatum Bud Appearance", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Dryas P2 - Appearance of white on first flower buds

# filter for only P2
drybud <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
drybud_anova <- aov(phase_DATE ~ obs, data = drybud)
summary(drybud_anova)

# boxplot to show distribution 2016-2019
(drybud_plot <- drybud %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First D. integrifolia Bud Appearance", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Dryas P3 - First open flower

# filter for only P3
dryopen <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P3") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
dryopen_anova <- aov(phase_DATE ~ obs, data = dryopen)
summary(dryopen_anova)

# boxplot to show distribution 2016-2019
(dryopen_plot <- dryopen %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First D. integrifolia open flower", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Dryas P4 - last petal shed

# filter for only P4
dryshed <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P4") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
dryshed_anova <- aov(phase_DATE ~ obs, data = dryshed)
summary(dryshed_anova)

# boxplot to show distribution 2016-2019
(dryshed_plot <- dryshed %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First D. integrifolia petal shed", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Dryas P5 - twisting of filament

# filter for only P5
drytwist <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P5") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
drytwist_anova <- aov(phase_DATE ~ obs, data = drytwist)
summary(drytwist_anova)

# boxplot to show distribution 2016-2019
(drytwist_plot <- drytwist %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "First D. integrifolia twisting of filament", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Salix P2 - first leaf bud burst

# filter for only P2
salbud <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
salbud_anova <- aov(phase_DATE ~ obs, data = salbud)
summary(salbud_anova)

# boxplot to show distribution 2016-2019
(salbud_plot <- salbud %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "S. arctica first leaf bud burst", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Salix P5 - first leaf yellowing

# filter for only P5
salsen1 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P5") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
salsen1_anova <- aov(phase_DATE ~ obs, data = salsen1)
summary(salsen1_anova)

# boxplot to show distribution 2016-2019
(salsen1_plot <- salsen1 %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "S. arctica first leaf turns yellow", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

# Salix P6 - all leaves turn yellow

# filter for only P6
salsen2 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P6") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# anova (h1 = sig dif between phenocam and transect obs, h0 = no sig dif between obs)
salsen2_anova <- aov(phase_DATE ~ obs, data = salsen2)
summary(salsen2_anova)

# boxplot to show distribution 2016-2019
(salsen2_plot <- salsen2 %>%
    ggplot() +
    aes(x = obs, y = phase_DATE, fill = obs) +
    geom_boxplot() +
    hrbrthemes::scale_fill_ipsum() +
    stat_compare_means(method = "anova") +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", title = "S. arctica last leaf turns yellow", fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none"))

#### PANEL OF ANOVA PLOTS ####

anova_pheno <- ggarrange(snow_free,eriobud_plot, drybud_plot,
                         dryopen_plot, dryshed_plot, drytwist_plot,
                         salbud_plot,salsen1_plot, salsen2_plot,
                         ncol = 3,nrow = 3)
anova_pheno

ggsave(anova_pheno, filename = "scripts/phenology_scripts/phenocams/anova_phenocam_box.png", height = 10, width = 12)

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
    xlim(125, 270) +
    theme_classic() +
    theme(legend.position = "none"))

#### all pheno phases chrono ###

(fulltimeline_chrono <- pheno %>%
  filter(!(Spp %in% c("SNOW",  "LUPARC", "PEDIC"))) %>%
  filter(!(phase_ID %in% c("S1", "S4","S5","S6","S2", "P1_before", "P2_before", "P3_before", "P4_before", 
                           "P5_before", "P6_before", "P7_before"))) %>%
  ggplot() +
  aes(x = Year, y = phase_DATE, colour = phase_ID) +
  geom_point(size = 1L) +
  stat_smooth(method = 'lm') +
  hrbrthemes::scale_color_ft() +
  theme_minimal() +
    labs(x = "Year", y = "DOY", fill = "Phenophase ID") +
    
  facet_wrap(vars(Spp)))
