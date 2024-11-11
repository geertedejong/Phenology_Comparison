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

#### LOAD FULL PHENOLOGY DATA DATA ####
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

#### Figure 1 comparative boxplots ####

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
                            "E.Vag.Snow Free" = "ERIVAG_P1",
                            "E.Vag. First Flower Bud" = "ERIVAG_P2",
                            "E.Vag. First Pollen visible" = "ERIVAG_P3",
                            "Dry.Int. Snow Free" = "DRYINT_P1",
                            "Dry.Int. First White on Bud" = "DRYINT_P2",
                            "Dry.Int. First Open Flower" = "DRYINT_P3",
                            "Dry.Int. Last Petal Shed" = "DRYINT_P4",
                            "Dry.Int. First Twisting of Filaments" = "DRYINT_P5",
                            "Dry.Int. All leaves dead" = "DRYINT_P6",
                            "Sal.Arc. Snow Free" = "SALARC_P1",
                            "Sal.Arc. First Leaf Bud Burst " = "SALARC_P2",
                            "Sal.Arc. First Pollen" = "SALARC_P3",
                            "Sal.Arc. Onset Seed Dispersal" = "SALARC_P4",
                            "Sal.Arc. First Yellow Leaf" = "SALARC_P5",
                            "Sal.Arc. Last Green Leaf" = "SALARC_P6",
                            "Sal.Arc. All Leaves Dead" = "SALARC_P7")


plot_boxplot <- function(data, x_var, y_var, x_limits = NULL, title = "Phenology stuff", xlabel = "DOY (2016 - 2019)", ylabel = "Phenophase") {
  data %>%
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }}, fill = interaction(Spp, obs))) +
    geom_boxplot() +
    #scale_fill_brewer(palette = "Set2") + # Use a qualitative palette that adapts to number of groups
    scale_fill_manual(values = c("lightgreen","darkgreen","yellow","orange","pink","purple","blue"),
                      breaks = c("SALARC.transect","SALARC.phenocam","DRYINT.transect","DRYINT.phenocam"  ,"ERIVAG.transect","ERIVAG.phenocam","SNOW.phenocam"),
                      labels = c("Sal.Arc. transect","Sal.Arc. phenocam","Dry.Int. transect","Dry.Int. phenocam","Eri.Vag. transect","Eri.Vag phenocam","Snow phenocam"))+
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

#### combination plot of cams, obs and NDVI ####
(comb_plot <- ggplot()+
   #geom_point(data=s2_ndvisf, aes(x=doi, y=NDVI_20m, color=factor(year)), alpha=0.3, size=1)+
   geom_smooth(data=s2_ndvisf,aes(x=doi, y=NDVI_20m, color=factor(year))) +
   hrbrthemes::scale_color_ipsum() +
   geom_vline(xintercept=cam_sf16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_senescence16, linetype='dashed',color='orange',size=1)+
   geom_vline(xintercept=cam_sf17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_senescence17, linetype='dashed',color='green',size=1)+
   geom_vline(xintercept=cam_sf18, linetype='dashed',color='purple',size=1)+
   geom_vline(xintercept=cam_senescence18, linetype='dashed',color='purple',size=1)+
   geom_vline(xintercept=cam_sf19, linetype='dashed',color='lightblue',size=1)+
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
   theme(legend.position = "right")
)

ggsave(comb_plot, filename = "figures/cam_obs_ndvi_greencurv.png", height = 10, width = 12)

#### Figure 3 big plot combining phenocams, transect and sat data ####

library(ggplot2)
library(dplyr)
library(hrbrthemes)

library(patchwork)

# Function to create the phenophase boxplot
plot_phenophase_boxplot <- function(data, x_var, y_var, x_limits = NULL, 
                                    title = "Phenology stuff", xlabel = "DOY (2016 - 2019)", ylabel = "Phenophase") {
  p <- data %>%
    ggplot() +
    geom_boxplot(aes(x = {{ x_var }}, y = {{ y_var }}, fill = interaction(Spp, obs))) +
    scale_fill_manual(values = c("lightgreen","darkgreen","yellow","orange","pink","purple","blue"),
                      breaks = c("SALARC.transect","SALARC.phenocam","DRYINT.transect","DRYINT.phenocam",
                                 "ERIVAG.transect","ERIVAG.phenocam","SNOW.phenocam"),
                      labels = c("Sal.Arc. transect","Sal.Arc. phenocam","Dry.Int. transect",
                                 "Dry.Int. phenocam","Eri.Vag. transect","Eri.Vag phenocam",
                                 "Snow phenocam")) +
    labs(x = xlabel, y = ylabel, title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "right") +
    coord_cartesian(xlim = x_limits)  # Ensure x-axis matches NDVI plot
  
  return(p)
}

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
phenophase_plot <- plot_phenophase_boxplot(overall2, phase_DATE, plot, x_limits = x_limits, title = "Phenophase Chronology")
ndvi_plot <- plot_ndvi_line(s2_ndvisf, x_limits = x_limits)

# Combine plots vertically using patchwork
combined_plot <- phenophase_plot / ndvi_plot

# Display the combined plot
combined_plot
