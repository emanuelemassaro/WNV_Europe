# APES: Land use and MBDs#######################################################
################################################################################
# The goal of this analysis is to explore potential associations between land
# cover indicators and the prevalence and incidence of MBDs at NUTS3 level in the
# whole larger European territory. 
# 
# To do this: 
# 1) First we considered averaged land use values (2006-19) and the overall 
# presence/absence of MBDs cases, per NUTS3 level both generally as well as by diseases. 
# This descriptive analysis was presented with box plots
#
# 2) Secondly, we considered the average land use (2006-19) and the total prevalence
# incidence of MBDs per NUTS3 level both generally as well as by diseases.
################################################################################

options(scipen=999) # Disable scientific notation

# Load & install the necessary packages
if(!require("pacman")) install.packages("pacman")

pacman::p_load(pacman, tidyverse, giscoR, sf, hrbrthemes, ggpubr, lubridate, expss)

#lc<- 
#  list.files(pattern = "*.csv") %>% 
#  map_df(~read_csv(.))
#lc <- lc[, -c(1,2)]
#saveRDS(lc, file = "lc_base.rds")

# LOAD THE DATA ################################################################
## Land cover
lc <- readRDS(file = "lc_base.rds")

## NUTS
### Codes
nuts_code <- read.csv("nuts_code.csv")

### Geometry
nuts <- giscoR::gisco_get_nuts(
  year = 2016,
  nuts_level = 3, 
  resolution = "03") %>%
  rename(NUTS_CODE = NUTS_ID) %>%
  filter(!grepl("FRY",NUTS_CODE))
  

## MBDs
### importing function
imp <- function(x) {
  tmp <- read_csv(x) %>%
    #filter(Imported == "N", LaboratoryResult == "CONF") %>%
    mutate(NUTS_CODE= ifelse(PlaceOfNotification != "UNK", PlaceOfNotification, PlaceOfInfection),
           Years = as.numeric(DateUsedForStatisticsYear)) %>%
    group_by(Years,
             NUTS_CODE,
             RecordType, Imported) %>%
    summarise(Total_cases = sum(NumberOfCases)) %>%
    as_tibble()
}

chik <- imp("CHIK_Case (1).csv")
deng <- imp("DENGE_Case (1).csv")
wnv <- imp("WNF_Case.csv")
yelf <- imp("YELF_Case (1).csv")
zik <- imp("ZIK_Case (1).csv")

# DATA PASS ####################################################################
## Create a common MBDs dataset and the following groups
mbds <- rbind(chik, deng, wnv, yelf, zik) %>% as_tibble()
mbds2 <- mbds %>% mutate(Cases = ifelse(is.na(Total_cases), 0, 1))

mbds.mean <- mbds %>% group_by(NUTS_CODE, Years) %>% summarise(Total_cases = sum(Total_cases))
mbds.mean2 <- mbds %>% group_by(NUTS_CODE, RecordType) %>% summarise(Total_cases = sum(Total_cases))

## Join different datasets
### Land cover and nuts code
nuts.full <- left_join(nuts, nuts_code, by = "NUTS_CODE") %>% as_tibble()
lc_nuts <- left_join(nuts.full, lc, by = "idx") %>% as_tibble()
colnames(lc_nuts) <- c("LEVL_CODE", "NUTS_CODE", "URBN_TYPE", "CNTR_CODE", "NAME_LATN", "NUTS_NAME", "MOUNT_TYPE", "COAST_TYPE", "FID", "geo", "idx", "geometry", "Crop_rain", "Herb_cover", "Tree_cover","Crop_irr", "Mosaic_crop", "Mosaic_natural", "Broadleaved_evergreen_15", "Broadleaved_decid_15", "Broadleaved_decid_40", 
                    "Broadleaved_decid_1540", "Needleleaved_evergreen_15", "Needleleaved_evergreen_40", "Needleleaved_evergreen_1540", "Needleleaved_decid_15", 
                    "Needleleaved_decid_40", "Needleleaved_decid_1540", "Mixleaves", "Mosaic_tree","Mosaic_herb", "Shrubland", "Shrubland_evergreen", "Shrubland_decid","Grassland", 
                    "Moss", "Sparse_vegetation", "Sparse_tree", "Sparse_shrub", "Sparse_herb", "Treecover_water", "Treecover_saltwater", "Shrub_water", "Urban_areas", "Bare_areas", 
                    "Bare_areas_cons", "Bare_areas_uncons", "Water", "Snow", "Years") 

lc_nuts.mean <- lc_nuts %>% 
  group_by(NUTS_CODE, Years) %>%
  summarise(Crop_rain = mean(Crop_rain),
                    Herb_cover = mean(Herb_cover),
                    Tree_cover = mean(Tree_cover),
                    Crop_irr = mean(Crop_irr),
                    Mosaic_crop = mean(Mosaic_crop),
                    Mosaic_natural = mean(Mosaic_natural),
                    Broadleaved_evergreen_15 = mean(Broadleaved_evergreen_15),
                    Broadleaved_decid_15 = mean(Broadleaved_decid_15),
                    Broadleaved_decid_40 = mean(Broadleaved_decid_40),
                    Broadleaved_decid_1540 = mean(Broadleaved_decid_1540),
                    Needleaved_evergreen_15 = mean(Needleleaved_evergreen_15),
                    Needleaved_evergreen_40 = mean(Needleleaved_evergreen_40),
                    Needleaved_evergreen_1540 = mean(Needleleaved_evergreen_1540),
                    Needleaved_decid_15 = mean(Needleleaved_decid_15),
                    Needleaved_decid_40 = mean(Needleleaved_decid_40),
                    Needleaved_decid_1540 = mean(Needleleaved_decid_1540),
                    Mixleaves = mean(Mixleaves),
                    Mosaic_tree = mean(Mosaic_tree),
                    Mosaic_herb = mean(Mosaic_herb),
                    Shrubland = mean(Shrubland),
                    Shrubland_evergreen = mean(Shrubland_evergreen),
                    Shrubland_decid = mean(Shrubland_decid),
                    Grassland = mean(Grassland),
                    Moss = mean(Moss),
                    Sparse_vegetation = mean(Sparse_vegetation),
                    Sparse_tree = mean(Sparse_tree),
                    Sparse_shrub = mean(Sparse_shrub),
                    Sparse_herb = mean(Sparse_herb),
                    Treecover_water = mean(Treecover_water),
                    Treecover_saltwater = mean(Treecover_saltwater),
                    Shrub_water = mean(Shrub_water),
                    Urban_areas = mean(Urban_areas),
                    Bare_areas = mean(Bare_areas),
                    Bare_areas_cons = mean(Bare_areas_cons),
                    Bare_areas_uncons = mean(Bare_areas_uncons),
                    Water = mean(Water),
                    Snow = mean(Snow)) %>%
  as_tibble()

## Land cover and NUTS geometry
#tmp <- left_join(nuts, lc_nuts.mean, by = "NUTS_CODE") %>% as_tibble()

## Land cover + geometry and MBDs
test <- left_join(lc_nuts, mbds.sum, by = c("NUTS_CODE", "Years")) %>% 
  mutate(cases = ifelse(is.na(Total_cases), 0, 1)) %>% # Add a categorical variable Cases ("Yes"/"No")
  as_tibble()

#test2 <- left_join(lc_nuts, mbds.sum2, by = "NUTS_CODE") %>% 
#  mutate(cases = ifelse(is.na(Total_cases), 0, 1)) %>% # Add a categorical variable Cases ("Yes"/"No")
#  as_tibble()

# LABEL THE FINAL DATASET ######################################################
colnames(test) <- c("LEVL_CODE", "NUTS_CODE", "URBN_TYPE", "CNTR_CODE", "NAME_LATN", "NUTS_NAME", "MOUNT_TYPE", "COAST_TYPE", "FID", "geo", "idx", "geometry", "Crop_rain", 
                    "Herb_cover", "Tree_cover","Crop_irr", "Mosaic_crop", "Mosaic_natural", "Broadleaved_evergreen_15", "Broadleaved_decid_15", "Broadleaved_decid_40", 
                    "Broadleaved_decid_1540", "Needleleaved_evergreen_15", "Needleleaved_evergreen_40", "Needleleaved_evergreen_1540", "Needleleaved_decid_15", 
                    "Needleleaved_decid_40", "Needleleaved_decid_1540", "Mixleaves", "Mosaic_tree","Mosaic_herb", "Shrubland", "Shrubland_evergreen", "Shrubland_decid","Grassland", 
                    "Moss", "Sparse_vegetation", "Sparse_tree", "Sparse_shrub", "Sparse_herb", "Treecover_water", "Treecover_saltwater", "Shrub_water", "Urban_areas", "Bare_areas", 
                    "Bare_areas_cons", "Bare_areas_uncons", "Water", "Snow", "Years", "Total_cases", "cases") 

#colnames(test2) <- c("LEVL_CODE", "NUTS_CODE", "URBN_TYPE", "CNTR_CODE", "NAME_LATN", "NUTS_NAME", "MOUNT_TYPE", "COAST_TYPE", "FID", "geo", "geometry", "Crop_rain", 
#                    "Herb_cover", "Tree_cover","Crop_irr", "Mosaic_crop", "Mosaic_natural", "Broadleaved_evergreen_15", "Broadleaved_decid_15", "Broadleaved_decid_40", 
#                    "Broadleaved_decid_1540", "Needleleaved_evergreen_15", "Needleleaved_evergreen_40", "Needleleaved_evergreen_1540", "Needleleaved_decid_15", 
#                    "Needleleaved_decid_40", "Needleleaved_decid_1540", "Mixleaves", "Mosaic_tree","Mosaic_herb", "Shrubland", "Shrubland_evergreen", "Shrubland_decid","Grassland", 
#                    "Moss", "Sparse_vegetation", "Sparse_tree", "Sparse_shrub", "Sparse_herb", "Treecover_water", "Treecover_saltwater", "Shrub_water", "Urban_areas", "Bare_areas", 
#                    "Bare_areas_cons", "Bare_areas_uncons", "Water", "Snow", "RecordType", "Total_cases") 

test <- apply_labels(test,
                     LEVL_CODE = "LEVL_CODE",
                     NUTS_CODE = "NUTS_CODE",
                     URBN_TYPE = "URBN_TYPE",
                     CNTR_CODE = "CNTR_CODE",
                     NAME_LATN = "NAME_LATN",
                     NUTS_NAME = "NUTS_NAME",
                     MOUNT_TYPE = "MOUNT_TYPE",
                     COAST_TYPE = "COAST_TYPE",
                     FID = "FID",
                     geo = "geo",
                     idx = "idx",
                     geometry = "geometry",
                     Crop_rain = "Cropland, rainfed",
                     Herb_cover = "Herbaceous cover",
                     Tree_cover = "Tree or shrub cover",
                     Crop_irr = "Cropland, irrigated or post‐flooding",
                     Mosaic_crop = "Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)",
                     Mosaic_natural = "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)",
                     Broadleaved_evergreen_15 = "Tree cover, broadleaved, evergreen, closed to open (>15%)",
                     Broadleaved_decid_15 = "Tree cover, broadleaved, deciduous, closed to open (>15%)",
                     Broadleaved_decid_40 = "Tree cover, broadleaved, deciduous, closed (>40%)",
                     Broadleaved_decid_1540 = "Tree cover, broadleaved, deciduous, open (15‐40%)",
                     Needleleaved_evergreen_15 = "Tree cover, needleleaved, evergreen, closed to open (>15%)",
                     Needleleaved_evergreen_40 = "Tree cover, needleleaved, evergreen, closed (>40%)",
                     Needleleaved_evergreen_1540 = "Tree cover, needleleaved, evergreen, open (15‐40%)",
                     Needleleaved_decid_15 = "Tree cover, needleleaved, deciduous, closed to open (>15%)",
                     Needleleaved_decid_40 = "Tree cover, needleleaved, deciduous, closed (>40%)",
                     Needleleaved_decid_1540 = "Tree cover, needleleaved, deciduous, open (15‐40%)",
                     Mixleaves = "Tree cover, mixed leaf type (broadleaved and needleleaved)",
                     Mosaic_tree = "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
                     Mosaic_herb = "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)",
                     Shrubland = "Shrubland",
                     Shrubland_evergreen = "Evergreen shrubland",
                     Shrubland_decid = "Deciduous shrubland",
                     Grassland = "Grassland",
                     Moss = "Lichens and mosses",
                     Sparse_vegetation = "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)",
                     Sparse_tree = "Sparse tree (<15%)",
                     Sparse_shrub = "Sparse shrub (<15%)",
                     Sparse_herb = "Sparse herbaceous cover (<15%)",
                     Treecover_water = "Tree cover, flooded, fresh or brakish water",
                     Treecover_saltwater = "Tree cover, flooded, saline water",
                     Shrub_water = "Shrub or herbaceous cover, flooded, fresh/saline/brakish water",
                     Urban_areas = "Urban areas",
                     Bare_areas = "Bare areas",
                     Bare_areas_cons = "Consolidated bare areas", 
                     Bare_areas_uncons = "Unconsolidated bare areas",
                     Water = "Water bodies",
                     Snow = "Permanent snow and ice",
                     Years = "Years", 
                     Total_cases = "Total cases",
                     cases = "cases yes/no") 

test2 <- apply_labels(test2,
                     LEVL_CODE = "LEVL_CODE",
                     NUTS_CODE = "NUTS_CODE",
                     URBN_TYPE = "URBN_TYPE",
                     CNTR_CODE = "CNTR_CODE",
                     NAME_LATN = "NAME_LATN",
                     NUTS_NAME = "NUTS_NAME",
                     MOUNT_TYPE = "MOUNT_TYPE",
                     COAST_TYPE = "COAST_TYPE",
                     FID = "FID",
                     geo = "geo",
                     geometry = "geometry",
                     Crop_rain = "Cropland, rainfed",
                     Herb_cover = "Herbaceous cover",
                     Tree_cover = "Tree or shrub cover",
                     Crop_irr = "Cropland, irrigated or post‐flooding",
                     Mosaic_crop = "Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)",
                     Mosaic_natural = "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)",
                     Broadleaved_evergreen_15 = "Tree cover, broadleaved, evergreen, closed to open (>15%)",
                     Broadleaved_decid_15 = "Tree cover, broadleaved, deciduous, closed to open (>15%)",
                     Broadleaved_decid_40 = "Tree cover, broadleaved, deciduous, closed (>40%)",
                     Broadleaved_decid_1540 = "Tree cover, broadleaved, deciduous, open (15‐40%)",
                     Needleleaved_evergreen_15 = "Tree cover, needleleaved, evergreen, closed to open (>15%)",
                     Needleleaved_evergreen_40 = "Tree cover, needleleaved, evergreen, closed (>40%)",
                     Needleleaved_evergreen_1540 = "Tree cover, needleleaved, evergreen, open (15‐40%)",
                     Needleleaved_decid_15 = "Tree cover, needleleaved, deciduous, closed to open (>15%)",
                     Needleleaved_decid_40 = "Tree cover, needleleaved, deciduous, closed (>40%)",
                     Needleleaved_decid_1540 = "Tree cover, needleleaved, deciduous, open (15‐40%)",
                     Mixleaves = "Tree cover, mixed leaf type (broadleaved and needleleaved)",
                     Mosaic_tree = "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
                     Mosaic_herb = "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)",
                     Shrubland = "Shrubland",
                     Shrubland_evergreen = "Evergreen shrubland",
                     Shrubland_decid = "Deciduous shrubland",
                     Grassland = "Grassland",
                     Moss = "Lichens and mosses",
                     Sparse_vegetation = "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)",
                     Sparse_tree = "Sparse tree (<15%)",
                     Sparse_shrub = "Sparse shrub (<15%)",
                     Sparse_herb = "Sparse herbaceous cover (<15%)",
                     Treecover_water = "Tree cover, flooded, fresh or brakish water",
                     Treecover_saltwater = "Tree cover, flooded, saline water",
                     Shrub_water = "Shrub or herbaceous cover, flooded, fresh/saline/brakish water",
                     Urban_areas = "Urban areas",
                     Bare_areas = "Bare areas",
                     Bare_areas_cons = "Consolidated bare areas", 
                     Bare_areas_uncons = "Unconsolidated bare areas",
                     Water = "Water bodies",
                     Snow = "Permanent snow and ice",
                     RecordType = "RecordType", 
                     Total_cases = "Total cases") 

# DATA ANALYSIS ################################################################
## Descriptive statistics
### Create function for plots
map.lc <- function(x, y) {
  tmp_string <- "Distribution ofin Europe"
  tmp_string2 <- gsub("^(.{15})(.*)$",         # Apply gsub
                      y,
                      tmp_string)
  
  ggplot(test)+
    geom_sf(aes(fill=log(x), geometry=geometry))+
    scale_fill_viridis_c(direction = -1,
                         na.value="white") +
    theme_ipsum()+
    theme(
      legend.position="right",
      plot.title = element_text(size=12))+
    labs(title = tmp_string2,
         subtile = "Land cover indicator is displayed aggregated by NUTS3-level, averaged over the period 2006-19",
         caption = "Source: ESA CCI Land Cover")
}

box.lc <- function(zx){
ggplot(test) +
  geom_boxplot(aes(x=cases, y=x))+
  theme_ipsum() +
  theme(
    legend.position="right",
    plot.title = element_text(size=12))+
  labs(title = tmp_string2,
       subtile = "Cases are displayed aggregated by NUTS3-level, and over the whole period of specific surveillance of the disease",
       caption = "Source: ECDC Dataset")
}

### Geographical mapping of the land cover indicators
map.lc(test$Crop_rain, "Cropland rainfed")
map.lc(test$Herb_cover, "Herb Cover")
map.lc(test$Tree_cover, "Tree or shrub cover")
map.lc(test$Crop_irr, "Cropland, irrigated or post‐flooding")
map.lc(test$Mosaic_crop, "Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)")
map.lc(test$Mosaic_natural, "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)")
map.lc(test$Broadleaved_evergreen_15, "Tree cover, broadleaved, evergreen, closed to open (>15%)")
map.lc(test$Broadleaved_decid_15, "Tree cover, broadleaved, deciduous, closed to open (>15%)")
map.lc(test$Broadleaved_decid_40, "Tree cover, broadleaved, deciduous, closed (>40%)")
map.lc(test$Broadleaved_decid_1540, "Tree cover, broadleaved, deciduous, open (15‐40%)")
map.lc(test$Needleleaved_evergreen_15, "Tree cover, needleleaved, evergreen, closed to open (>15%)")
map.lc(test$Needleleaved_evergreen_40, "Tree cover, needleleaved, evergreen, closed (>40%)")
map.lc(test$Needleleaved_evergreen_1540, "Tree cover, needleleaved, evergreen, open (15‐40%)")
map.lc(test$Needleleaved_decid_15, "Tree cover, needleleaved, deciduous, closed to open (>15%)")
map.lc(test$Needleleaved_decid_40, "Tree cover, needleleaved, deciduous, closed (>40%)")
map.lc(test$Needleleaved_decid_1540, "Tree cover, needleleaved, deciduous, open (15‐40%)")
map.lc(test$Mixleaves, "Tree cover, mixed leaf type (broadleaved and needleleaved)")
map.lc(test$Mosaic_tree, "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)")
map.lc(test$Mosaic_herb, "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)")
map.lc(test$Shrubland, "Shrubland")
map.lc(test$Shrubland_evergreen, "Evergreen shrubland")
map.lc(test$Shrubland_decid, "Deciduous shrubland")
map.lc(test$Grassland, "Grassland")
map.lc(test$Moss, "Lichens and mosses")
map.lc(test$Sparse_vegetation, "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)")
map.lc(test$Sparse_tree, "Sparse tree (<15%)")
map.lc(test$Sparse_shrub, "Sparse shrub (<15%)")
map.lc(test$Sparse_herb, "Sparse herbaceous cover (<15%)")
map.lc(test$Treecover_water, "Tree cover, flooded, fresh or brakish water")
map.lc(test$Treecover_saltwater, "Tree cover, flooded, saline water")
map.lc(test$Shrub_water, "Shrub or herbaceous cover, flooded, fresh/saline/brakish water")
map.lc(test$Urban_areas, "Urban areas")
map.lc(test$Bare_areas, "Bare areas")
map.lc(test$Bare_areas_cons, "Consolidated bare areas") 
map.lc(test$Bare_areas_uncons, "Unconsolidated bare areas")
map.lc(test$Water, "Water bodies")
map.lc(test$Snow, "Permanent snow and ice")
map.lc(test$RecordType, "RecordType")
map.lc(test$tot_cases, "tot_cases")
map.lc(test$cases, "cases yes/no") 


### Boxplots
#### BY CASES YES/NO
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Crop_rain))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Herb_cover))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Tree_cover))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Crop_irr))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Mosaic_crop))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Mosaic_natural))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Broadleaved_evergreen_15))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Broadleaved_decid_15))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Broadleaved_decid_40))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Broadleaved_decid_1540))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_evergreen_15))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_evergreen_40))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_evergreen_1540))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_decid_15))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_decid_40))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Needleleaved_decid_1540))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Mixleaves))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Mosaic_tree))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Mosaic_herb))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Shrubland))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Shrubland_evergreen))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Shrubland_decid))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Grassland))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Moss))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Sparse_vegetation))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Sparse_tree))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Sparse_shrub))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Sparse_herb))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Treecover_water))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Treecover_saltwater))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Shrub_water))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Urban_areas))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Bare_areas))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Bare_areas_cons))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Bare_areas_uncons))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Water))) +theme_ipsum()
ggplot(test) +
  geom_boxplot(aes(x=cases, y=log(Snow))) +theme_ipsum()


#### BY RECORD TYPE
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Crop_rain))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Herb_cover))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Tree_cover))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Crop_irr))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Mosaic_crop))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Mosaic_natural))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Broadleaved_evergreen_15))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Broadleaved_decid_15))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Broadleaved_decid_40))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Broadleaved_decid_1540))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_evergreen_15))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_evergreen_40))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_evergreen_1540))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_decid_15))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_decid_40))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Needleleaved_decid_1540))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Mixleaves))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Mosaic_tree))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Mosaic_herb))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Shrubland))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Shrubland_evergreen))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Shrubland_decid))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Grassland))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Moss))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Sparse_vegetation))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Sparse_tree))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Sparse_shrub))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Sparse_herb))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Treecover_water))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Treecover_saltwater))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Shrub_water))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Urban_areas))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Bare_areas))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Bare_areas_cons))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Bare_areas_uncons))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Water))) +theme_ipsum()
ggplot(test2) +
  geom_boxplot(aes(x=RecordType, y=log(Snow))) +theme_ipsum()

## Testing the differences in mean land cover parameters between NUTS3 regions with
## and without cases
t.test(Crop_rain ~ cases, data = test, var.equal = TRUE)
t.test(Herb_cover ~ cases, data = test, var.equal = TRUE) #p-value = 0.0318
t.test(Tree_cover ~ cases, data = test, var.equal = TRUE) 
t.test(Crop_irr ~ cases, data = test, var.equal = TRUE)
t.test(Mosaic_crop ~ cases, data = test, var.equal = TRUE)
t.test(Mosaic_natural ~ cases, data = test, var.equal = TRUE) #p-value = 0.000001725
t.test(Broadleaved_evergreen_15 ~ cases, data = test, var.equal = TRUE)
t.test(Broadleaved_decid_15 ~ cases, data = test, var.equal = TRUE) #p-value = 0.000174
t.test(Broadleaved_decid_40 ~ cases, data = test, var.equal = TRUE)
t.test(Broadleaved_decid_1540 ~ cases, data = test, var.equal = TRUE)
t.test(Needleleaved_evergreen_15 ~ cases, data = test, var.equal = TRUE) #p-value =  0.00006788
t.test(Needleleaved_evergreen_40 ~ cases, data = test, var.equal = TRUE)
t.test(Needleleaved_evergreen_1540 ~ cases, data = test, var.equal = TRUE)
t.test(Needleleaved_decid_15 ~ cases, data = test, var.equal = TRUE)
t.test(Needleleaved_decid_40 ~ cases, data = test, var.equal = TRUE)
t.test(Needleleaved_decid_1540 ~ cases, data = test, var.equal = TRUE)
t.test(Mixleaves ~ cases, data = test, var.equal = TRUE)
t.test(Mosaic_tree ~ cases, data = test, var.equal = TRUE) #p-value = 0.0000006264
t.test(Mosaic_herb ~ cases, data = test, var.equal = TRUE) #p-value = 0.00004644
t.test(Shrubland ~ cases, data = test, var.equal = TRUE)
t.test(Shrubland_evergreen ~ cases, data = test, var.equal = TRUE)
t.test(Shrubland_decid ~ cases, data = test, var.equal = TRUE) #p-value = 0.01836
t.test(Grassland ~ cases, data = test, var.equal = TRUE)
t.test(Moss ~ cases, data = test, var.equal = TRUE)
t.test(Sparse_vegetation ~ cases, data = test, var.equal = TRUE)
t.test(Sparse_tree ~ cases, data = test, var.equal = TRUE)
t.test(Sparse_shrub ~ cases, data = test, var.equal = TRUE) #p-value = 0.0265
t.test(Sparse_herb ~ cases, data = test, var.equal = TRUE) #p-value = 0.03019
t.test(Treecover_water ~ cases, data = test, var.equal = TRUE) #p-value = 0.0005618
t.test(Treecover_saltwater ~ cases, data = test, var.equal = TRUE)
t.test(Shrub_water ~ cases, data = test, var.equal = TRUE) #p-value = 0.006737
t.test(Urban_areas ~ cases, data = test, var.equal = TRUE) #p-value = 2.221e-10
t.test(Bare_areas ~ cases, data = test, var.equal = TRUE)
t.test(Bare_areas_cons ~ cases, data = test, var.equal = TRUE) #p-value = 0.0129
t.test(Bare_areas_uncons ~ cases, data = test, var.equal = TRUE)
t.test(Water ~ cases, data = test, var.equal = TRUE) #p-value = 0.0003265
t.test(Snow ~ cases, data = test, var.equal = TRUE)


## Testing same differences by diseases group
#test2 %<>% mutate(Disease = ifelse(is.na(RecordType), "No disease", RecordType)) %>% as_tibble()
addmargins(table(test2$RecordType, useNA = "always"))

anova <- aov(Crop_rain ~ RecordType, data = test2) #<0.0000000000000002
summary(anova)
n(anova)
anova <- aov(Herb_cover~ Disease, data = test2)
summary(anova)
anova <- aov(Tree_cover ~ Disease, data = test2) #0.0318
summary(anova)
anova <- aov(Crop_irr ~ Disease, data = test2) # 0.000000000462
summary(anova)
anova <- aov(Mosaic_crop ~ Disease, data = test2)
summary(anova)
anova <- aov(Mosaic_natural ~ Disease, data = test2) #0.00000000000237
summary(anova)
anova <- aov(Broadleaved_evergreen_15 ~ Disease, data = test2)
summary(anova)
anova <- aov(Broadleaved_decid_15 ~ Disease, data = test2) #<0.0000000000000002
summary(anova)
anova <- aov(Broadleaved_decid_40 ~ Disease, data = test2) #0.0000574
summary(anova)
anova <- aov(Broadleaved_decid_1540 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_evergreen_15 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_evergreen_40 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_evergreen_1540 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_decid_15 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_decid_40 ~ Disease, data = test2)
summary(anova)
anova <- aov(Needleleaved_decid_1540 ~ Disease, data = test2)
summary(anova)
anova <- aov(Mixleaves ~ Disease, data = test2)
summary(anova)
anova <- aov(Mosaic_tree ~ Disease, data = test2) #0.0000000055
summary(anova)
anova <- aov(Mosaic_herb ~ Disease, data = test2) #0.000018
summary(anova)
anova <- aov(Shrubland ~ Disease, data = test2) #0.000000623
summary(anova)
anova <- aov(Shrubland_evergreen ~ Disease, data = test2)
summary(anova)
anova <- aov(Shrubland_decid ~ Disease, data = test2)
summary(anova)
anova <- aov(Grassland ~ Disease, data = test2) #0.0035
summary(anova)
anova <- aov(Moss ~ Disease, data = test2)
summary(anova)
anova <- aov(Sparse_vegetation ~ Disease, data = test2)
summary(anova)
anova <- aov(Sparse_tree ~ Disease, data = test2)
summary(anova)
anova <- aov(Sparse_shrub ~ Disease, data = test2)
summary(anova)
anova <- aov(Sparse_herb ~ Disease, data = test2)
summary(anova)
anova <- aov(Treecover_water ~ Disease, data = test2) #0.00081
summary(anova)
anova <- aov(Treecover_saltwater ~ Disease, data = test2)
summary(anova)
anova <- aov(Shrub_water ~ Disease, data = test2) #0.0162
summary(anova)
anova <- aov(Urban_areas ~ Disease, data = test2) #<0.0000000000000002
summary(anova)
anova <- aov(Bare_areas ~ Disease, data = test2)
summary(anova)
anova <- aov(Bare_areas_cons ~ Disease, data = test2) #0.0264
summary(anova)
anova <- aov(Bare_areas_uncons ~ Disease, data = test2) #0.0276
summary(anova)
anova <- aov(Water ~ Disease, data = test2, na.rm = T) #0.0000854
summary(anova)
anova <- aov(Snow ~ Disease, data = test2, na.rm = T)
summary(anova)

