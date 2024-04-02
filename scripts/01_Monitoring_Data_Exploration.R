### 1. Phenology Monitoring Exploration and Preliminary Figures ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ####
### Date: 9th November 2020 ###

#### LOAD PACKAGES  #####

library(dplyr)
library(readr)
library(tidyverse) 
library(esquisse)

#### LOAD UP TO DATE TRANSECT DATA ####
transect <- read_csv(file = "data/qiki_phen_with_before_2019.txt")
str(transect)

# Duplicate plot ID and call it 'sex'
transect$sex = transect$Plot.ID

# keep sex marker (but not plot ID number) in sex column
transect$sex <- sub("^([[:alpha:]]*).*", "\\1", transect$sex)

# keep plot ID number (but not sex) in Plot.ID column
transect$Plot.ID <- parse_number(transect$Plot.ID)

#### LOAD & CLEAN PHENOCAM DATA ####
snow <- read_csv(file = "data/cam_SNOW.csv")
dryint <- read_csv(file = "data/cam_DRYINT.csv")
erivag <- read_csv(file = "data/cam_ERIVAG.csv")
salarc <- read_csv(file = "data/cam_SALARC.csv")
luparc <- read_csv(file = "data/cam_LUPARC.csv")
pedic <- read_csv(file = "data/cam_PEDIC.csv")

# the resulting data frames are not of the same size
# add columns containing NA to allow matching:

snow[, c("P2", "P2_Cert", "P3","P3_Cert","P4","P4_Cert","P5","P5_Cert","P6","P6_Cert",
         "Q1", "Q1_Cert", "Q2","Q2_Cert","Q3-M","Q3-M_Cert","Q3-F","Q3-F_Cert")] <- NA
dryint[, c("S1","S2","S3","S4","S5","S6","Q3-M","Q3-M_Cert","Q3-F","Q3-F_Cert")] <- NA
erivag[, c("P3","P3_Cert","P4","P4_Cert","P5","P5_Cert","P6","P6_Cert",
           "Q1", "Q1_Cert","Q3-M","Q3-M_Cert","Q3-F","Q3-F_Cert","S1","S2","S3","S4","S5","S6")] <- NA
salarc[, c("P3","P3_Cert","P4","P4_Cert","Q2","Q2_Cert",
           "Q1", "Q1_Cert","S1","S2","S3","S4","S5","S6")] <- NA
luparc[, c("P4","P4_Cert","P5","P5_Cert","P6","P6_Cert",
           "Q1", "Q1_Cert","Q3-M","Q3-M_Cert","Q3-F","Q3-F_Cert","S1","S2","S3","S4","S5","S6")] <- NA
pedic[, c("P4","P4_Cert","P5","P5_Cert","P6","P6_Cert",
          "Q1", "Q1_Cert","Q3-M","Q3-M_Cert","Q3-F","Q3-F_Cert","S1","S2","S3","S4","S5","S6")] <- NA

# Combine all into one  phenolocam data set
cam_phen <- rbind(snow,dryint,erivag,salarc,pedic,luparc)

# convert the three columns 'SPP' 'Pheno_ID' and 'year' stage into factors
cam_phen[,1:3] <- lapply(cam_phen[,1:3], factor)
transect[,1:3] <- lapply(transect[,1:3], factor)

str(cam_phen)

# remove all incomplete data entries that do not match the format DD/MM/YYYY
# use is.na() and grepl using a regular expression 
is.na(cam_phen$S1) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S1)
is.na(cam_phen$S2) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S2)
is.na(cam_phen$S3) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S3)
is.na(cam_phen$S4) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S4)
is.na(cam_phen$S5) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S5)
is.na(cam_phen$S6) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$S6)
is.na(cam_phen$P2) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$P2)
is.na(cam_phen$P3) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$P3)
is.na(cam_phen$P4) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$P4)
is.na(cam_phen$P5) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$P5)
is.na(cam_phen$P6) <- !grepl("\\d{2}/{1}\\d{2}/{1}\\d{4}", cam_phen$P6)

# convert the dates into day of year (DOY) by converting the date columns as a date 
cam_phen$S1 <- as.Date(cam_phen$S1, "%d/%m/%Y")
cam_phen$S2 <- as.Date(cam_phen$S2, "%d/%m/%Y")
cam_phen$S3 <- as.Date(cam_phen$S3, "%d/%m/%Y")
cam_phen$S4 <- as.Date(cam_phen$S4, "%d/%m/%Y")
cam_phen$S5 <- as.Date(cam_phen$S5, "%d/%m/%Y")
cam_phen$S6 <- as.Date(cam_phen$S6, "%d/%m/%Y")
cam_phen$P2 <- as.Date(cam_phen$P2, "%d/%m/%Y")
cam_phen$P3 <- as.Date(cam_phen$P3, "%d/%m/%Y")
cam_phen$P4 <- as.Date(cam_phen$P4, "%d/%m/%Y")
cam_phen$P5 <- as.Date(cam_phen$P5, "%d/%m/%Y")
cam_phen$P6 <- as.Date(cam_phen$P6, "%d/%m/%Y")

# turn dates into day of year (DOY)

cam_phen$S1 <- as.numeric(format(cam_phen$S1, "%j"))
cam_phen$S2 <- as.numeric(format(cam_phen$S2, "%j"))
cam_phen$S3 <- as.numeric(format(cam_phen$S3, "%j"))
cam_phen$S4 <- as.numeric(format(cam_phen$S4, "%j"))
cam_phen$S5 <- as.numeric(format(cam_phen$S5, "%j"))
cam_phen$S6 <- as.numeric(format(cam_phen$S6, "%j"))
cam_phen$P2 <- as.numeric(format(cam_phen$P2, "%j"))
cam_phen$P3 <- as.numeric(format(cam_phen$P3, "%j"))
cam_phen$P4 <- as.numeric(format(cam_phen$P4, "%j"))
cam_phen$P5 <- as.numeric(format(cam_phen$P5, "%j"))
cam_phen$P6 <- as.numeric(format(cam_phen$P6, "%j"))

# rename column SPP to Spp to match transect data
names(cam_phen)[names(cam_phen) == "SPP"] <- "Spp"
names(cam_phen)[names(cam_phen) == "Pheno_ID"] <- "Plot.ID"

#### MERGE TRANSECT & PHENOCAM DATA ####
# new column with observation type for each df
transect$obs <- 'transect'
cam_phen$obs <- 'phenocam'

# rename phenocam sites to match transects
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "1" = "Phenocam 1")) 
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "2" = "Phenocam 2")) 
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "3" = "Phenocam 3")) 
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "4" = "Phenocam 4")) 
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "5" = "Phenocam 5")) 
cam_phen <- cam_phen %>% 
  mutate(cam_phen, Plot.ID = fct_recode(Plot.ID, "6" = "Phenocam 6")) 
str(cam_phen)


# full join these tables using shared keys
pheno_full <- full_join(cam_phen, transect,copy=TRUE)
str(pheno_full)


#### LOAD PHENOCAM COORDS ####
pheno_coords <- read_csv(file = "data/phenocam_coordinates.csv")
pheno_coords <- pheno_coords %>% 
  mutate(Plot.ID = as.factor(Plot.ID))
pheno_full <- full_join(pheno_full, pheno_coords,copy=TRUE)
str(pheno_full)

# concatenate ID
pheno_full$ind.ID <- paste(pheno_full$Spp, pheno_full$Plot.ID, pheno_full$Year, sep="_")

# elongate the dataset so theres a phenophase column

pheno_long <- gather(pheno_full, phase_ID, phase_DATE,                           # in this order: data frame, key, value
                          c(S1, S2, S3, S4, S5, S6,P1,P2,P3,P4,P5,P6,P7,
                            P1_before,P2_before,P3_before,P4_before,P5_before,P6_before,P7_before))        

pheno_long <- gather(pheno_long, Q_ID, Q,                           # in this order: data frame, key, value
                     c(Q1, Q2, 'Q3-M', 'Q3-F'))        # we need to specify which columns to gather

pheno_long <- gather(pheno_long, cert_ID, cert,                           # in this order: data frame, key, value
                     c('P2_Cert','P3_Cert','P4_Cert',
                       'P5_Cert','P6_Cert','Q1_Cert', 'Q2_Cert', 'Q3-M_Cert', 'Q3-F_Cert'))        # we need to specify which columns to gather


# make sure all phenophase data is same type
pheno_long$phase_DATE <- as.numeric(pheno_long$phase_DATE)

# change year to a num
# save pheno_full dataset
#### Save full prop dataset as its own CSV ####
write.csv(pheno_long, file = "data/phenology_transect_cam.csv", row.names = FALSE)

####  GRAPH EACH SPECIES   ####

# Plot 1 - Snow Melt Date for all species 2001-2019 #
snowmelt_plot <- ggplot(transect) +
  aes(x = P3, fill = Spp) +
  geom_density(adjust = 1L, alpha = .7) +
  hrbrthemes::scale_fill_ipsum() +
  hrbrthemes::scale_color_ipsum() +
  labs(x = "Snow Melt Day of Year - Average 2001 - 2019", y = "Density") +
  theme_classic() +
  theme(legend.position = "bottom")

(snowmelt_plot <- snowmelt_plot + labs(fill = "Species"))


# Plot 2 - Bud burst dat for all species 2001-2019 #
(budburst_plot <- ggplot(pheno_full) +
    aes(x = P3, fill = Spp) +
    geom_density(adjust = 1L, alpha = .7) +
    hrbrthemes::scale_fill_ipsum() +
    hrbrthemes::scale_color_ipsum() +
    labs(x = "Bud Burst Day of Year - Average 2001 - 2019", y = "Density") +
    theme_classic() +
    theme(legend.position = "bottom"))

(budburst_plot <- budburst_plot + labs(fill = "Species"))


(snowmelt_v_budburst <- ggplot(transect) +
  aes(x = P1, y = P2, colour = Spp) +
  geom_point(size = 1L) +
  geom_smooth(span = 0.75) +
  hrbrthemes::scale_fill_ipsum() +
  hrbrthemes::scale_color_ipsum() +
  labs(x = "Date of snow melt (DOY)", y = "Date of bud burst (DOY)", color = "Species") +
  theme_minimal())

(snowmelt_chronology <- ggplot(transect) +
  aes(x = Year, y = P1) +
  geom_point(size = 1L, colour = "#0c4c8a") +
  geom_smooth(span = 0.75) +
  theme_minimal())

