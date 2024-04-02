### 3. Compare transect and phenocam data with satellite data ###
### Geerte de Jong, geerte.de.jong@bioenv.gu.se ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ###
### Date: 16th March 2021 ###

#### LOAD PACKAGES  #####

library(dplyr)
library(readr)
library(tidyverse) 
library(esquisse)
library(ggpubr)
library(gridExtra)


#### LOAD DATA ####
pheno <- read.csv(file = "data/phenology_transect_cam.csv")
s2    <- read.csv(file = "data/S2QHIphenocam.csv")
ndvi_m <- read.csv(file = "data/NDVI_modis.csv")
ndsi_m <- read.csv(file = "data/NDSI_modis.csv")

str(pheno)

#### Data wrangling copied from 02_Phenology_Obs_Comparison script written by Elise ####
# rename column s3 to p1 to match transect data
pheno <- pheno %>% 
  mutate(pheno, phase_ID = fct_recode(phase_ID, "P1" = "S3")) 

# filter for only P1
snowfree <- pheno %>% filter(phase_ID == "P1") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P2
eriobud <- pheno %>% filter(Spp %in% "ERIVAG") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P2
drybud <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P3
dryopen <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P3") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P4
dryshed <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P4") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P5
drytwist <- pheno %>% filter(Spp %in% "DRYINT") %>% 
  filter(phase_ID == "P5") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P2
salbud <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P2") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P5
salsen1 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P5") %>%
  filter(Year >= 2016L & Year <= 2019L) 

# filter for only P5
salsen2 <- pheno %>% filter(Spp %in% "SALARC") %>% 
  filter(phase_ID == "P6") %>%
  filter(Year >= 2016L & Year <= 2019L) 





#### some cleaning ####
s2_ndvisf<- subset(s2, NDVI_20m>0.2) #remove all NDVI values below o.2 to exclude negatives and snow
s2_ndsi <- subset(s2, NDSI_20m>0.4)

m_ndvisf<- subset(ndvi_m, NDVI>0.2) #remove all NDVI values below o.2 to exclude negatives and snow
m_ndsi <- subset(ndsi_m, NDSI>4000)
m_ndsi$NDSI <- (m_ndsi$NDSI)*0.0001

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
   theme(legend.position = "none"))

# average geom smooth and colors for camera locations, average for all years
(s2_plot <- s2_ndvisf %>%
    ggplot() +
    aes(x = doi, y = NDVI_20m) +
    geom_smooth() +
    geom_point(aes(color=year))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
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

#### Exploration plot of MODIS data ####
# NDVI
# average geom smooth and colors for camera locations, separate years
(m_plot <- m_ndvisf %>%
   ggplot() +
   aes(x = doy, y = NDVI) +
   geom_smooth() +
   geom_point(aes(color=name))+
   hrbrthemes::scale_fill_ipsum() +
   labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
   facet_grid(year~.) +
   theme_classic() +
   theme(legend.position = "none"))

# average geom smooth and colors for camera locations, average for all years
(m_plot <- m_ndvisf %>%
    ggplot() +
    aes(x = doy, y = NDVI) +
    geom_smooth() +
    geom_point(aes(color=year))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
    theme_classic() +
    theme(legend.position = "none"))

# NDSI
# average geom smooth and colors for camera locations, separate years
(m_plot <- m_ndsi %>%
    ggplot() +
    aes(x = doy, y = NDSI) +
    geom_smooth() +
    geom_point(aes(color=name))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDSI", x = "DOY (2016 - 2019)", fill = "year") +
    facet_grid(year~.) +
    theme_classic() +
    theme(legend.position = "none"))

# average geom smooth and colors for camera locations, average for all years
(m_plot <- m_ndsi %>%
    ggplot() +
    aes(x = doy, y = NDSI) +
    geom_smooth() +
    geom_point(aes(color=year))+
    hrbrthemes::scale_fill_ipsum() +
    labs(y = "NDSI", x = "DOY (2016 - 2019)", fill = "year") +
    theme_classic() +
    theme(legend.position = "none"))

#### combination plot S2 and MODIS ####
comb_db <- rbind(m_ndsi, m_ndvisf, s2_ndvisf, s2_ndsi)

# NB run 02_Phenology_Obs_Comparison script first to load data into environment
(cam_sf <- min(snowfree$phase_DATE, na.rm=TRUE))
(cam_greenup <- min(dryopen$phase_DATE, na.rm=TRUE))
(cam_senescence <- max(salsen2$phase_DATE, na.rm=TRUE))

#### combination plot of cams, obs and NDVI ####
(comb_plot <- ggplot()+
  geom_point(data=s2_ndvisf, aes(x=doi, y=NDVI_20m, color='coral1'), alpha=0.3, size=1 ,inherit.aes = FALSE)+
  geom_point(data=m_ndvisf, aes(x=doy, y=NDVI, color='deepskyblue'), alpha=0.3, size=1,inherit.aes = FALSE)+
  geom_smooth(data=s2_ndvisf,aes(x=doi, y=NDVI_20m, color='coral1'),inherit.aes = FALSE) +
  geom_smooth(data=m_ndvisf, aes(x=doy, y=NDVI, color='deepskyblue'),inherit.aes = FALSE) +
  hrbrthemes::scale_fill_ipsum() +
  geom_vline(xintercept=cam_sf, linetype='dashed')+
  geom_vline(xintercept=cam_greenup, linetype='dashed')+
  geom_vline(xintercept=cam_senescence, linetype='dashed')+
  xlim(100,300)+
  ylim(0.2,1)+
  labs(y = "NDVI", x = "DOY (2016 - 2019)", fill = "year") +
  theme_classic() +
  theme(legend.position = "none")
 )

ggsave(comb_plot, filename = "figures/cam_obs_ndvi_greencurv.png", height = 10, width = 12)

#...and NDSI
(comb_plot <- ggplot()+
    geom_point(data=s2_ndsi, aes(x=doi, y=NDSI_20m, color='coral1'), alpha=0.3, size=1 ,inherit.aes = FALSE)+
    geom_point(data=m_ndsi, aes(x=doy, y=NDSI, color='deepskyblue'), alpha=0.3, size=1,inherit.aes = FALSE)+
    geom_smooth(data=s2_ndsi,aes(x=doi, y=NDSI_20m, color='coral1'),inherit.aes = FALSE) +
    geom_smooth(data=m_ndsi, aes(x=doy, y=NDSI, color='deepskyblue'),inherit.aes = FALSE) +
    hrbrthemes::scale_fill_ipsum() +
    geom_vline(xintercept=cam_sf, linetype='dashed')+
    geom_vline(xintercept=cam_senescence, linetype='dashed')+
    ylim(0.4,1)+
    labs(y = "NDSI", x = "DOY (2016 - 2019)", fill = "year",inherit.aes = FALSE) +
    theme_classic() +
    theme(legend.position = "none")
)

ggsave(comb_plot, filename = "figures/cam_obs_ndsi_greencurv.png", height = 10, width = 12)
