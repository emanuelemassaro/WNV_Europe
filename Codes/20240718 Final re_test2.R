# LOAD THE LIBRARIES ###########################################################
library(corrplot)
library(gridExtra)
library(GWmodel)
library(MASS)
library(viridis)
library(spdep)
library(sf)
library(terra)
library(tidyverse)
library(xlsx)
library(hrbrthemes)

#
# SET WORKSPACE ###########################################################
options(prompt="R> ", digits=4, scipen=999) 
setwd("~/APES/Analyses/MBDs GWR")

#
# LOAD DATASET AND SELECT VARIABLES FOR MODEL##############################
corr.test <- readRDS("20240627_APES FINAL_metric.rds") %>%
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") %>%  
  dplyr::select(4, 6:42, 44) %>%
  group_by(geometry) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::select(-(geometry)) %>%
  na.omit() %>%
  as.matrix()

corr.matrix <- round(cor(corr.test), 2)
#xlsx::write.xlsx(as.data.frame(corr.matrix), "correlation_matrix.xlsx")
#
# Selected variables: Forest, Shrub, Urban, Crop, Water, Other, Max_WD_Winter, 
# Max_DD_Winter, Mean_T_Spring, Mean_P_Spring, Max_DD_Spring, Mean_T_Summer, 
# Mean_P_Autumn, Max_WD_Autumn, Max_DD_Autumn, gdp, pop_dens
corr.plot <- corrplot(corr.matrix, method = 'color', order = 'alphabet')
ggsave(corr.plot, filename = "Corrplot_APES.png", bg = "white", dpi=300)

# 
#
# LOAD AND PREP THE DATASET ####################################################
avg_var <- c("Forest", "Shrub", "Urban", "Crop", "Water", "Other", 
             "Mean_P_Winter", "Max_WD_Winter", "Max_DD_Winter", 
             "Mean_T_Spring", "Max_DD_Spring", 
             "Mean_T_Summer",
             "Max_WD_Autumn", "Max_DD_Autumn", 
             "gdp", "pop_dens")

fin.eu3 <- readRDS("20240627_APES FINAL_metric.rds") |> 
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") |>
   dplyr::select(incidence, Forest, Shrub, Urban, Crop, Water, Other, 
                 Mean_P_Winter, Max_WD_Winter, Max_DD_Winter, 
                 Mean_T_Spring, Max_DD_Spring, 
                 Mean_T_Summer, 
                 Max_WD_Autumn, Max_DD_Autumn, 
                 gdp, pop_dens, geometry, NUTS_NAME, NUTS_CODE) |>
  group_by(geometry, NUTS_NAME, NUTS_CODE) |>
  dplyr::summarise(incidence = sum(incidence, na.rm = TRUE), 
            across(all_of(avg_var), mean, na.rm = TRUE)) |>
  mutate(cent = sf::st_centroid(geometry)) |>
  na.omit() |>
  as_tibble()

# Standardizing the dataset
## Selecting variables to standardize
vars_to_scale <- fin.eu3 %>%
  st_drop_geometry() %>%  # Drop geometry column
  select_if(is.numeric) %>% # Select numeric columns
  colnames()

## Standardizing selected variables
scaled_data <- fin.eu3 %>%
  mutate(across(all_of(vars_to_scale), scale)) 

## Creating a weight matrix based on metric distance to allow for spatial-lag
## Changing from Euclydean to metric distance
### Creating an sf object
tmp <- sf::st_as_sf(scaled_data) # Create an sf object
scaled_data_meter <- na.omit(st_transform(tmp, crs = 32633)) # Transform it to metric
print(st_crs(scaled_data_meter)) # Check if is in metric

### Creating coordinates and max distance
scaled_spatial <- as(scaled_data_meter, "Spatial")
coords <- coordinates(scaled_spatial)
max_distance <- 150000

### Creeting a distance matrix
w <- dnearneigh(coords, 0, max_distance, longlat = FALSE)
w <- nb2listw(w, style = "B", zero.policy = TRUE)
W_matrix <- listw2mat(w)

### Moran's I
# List of variables you want to calculate Moran's I for
morans_incidence <- moran.test(fin.eu3$incidence, listw = w) %>% print()
morans_forest <- moran.test(fin.eu3$Forest, listw = w) %>% print()
morans_Shrub <- moran.test(fin.eu3$Shrub, listw = w) %>% print()
morans_Urban <- moran.test(fin.eu3$Urban, listw = w) %>% print()
morans_Crop <- moran.test(fin.eu3$Crop, listw = w) %>% print()
morans_Water <- moran.test(fin.eu3$Water, listw = w) %>% print()
morans_Other <- moran.test(fin.eu3$Other, listw = w) %>% print()
morans_Mean_P_Winter <- moran.test(fin.eu3$Mean_P_Winter, listw = w) %>% print()
morans_Max_WD_Winter <- moran.test(fin.eu3$Max_WD_Winter, listw = w) %>% print()
morans_Max_DD_Winter <- moran.test(fin.eu3$Max_DD_Winter, listw = w) %>% print()
morans_Mean_T_Spring <- moran.test(fin.eu3$Mean_T_Spring, listw = w) %>% print()
morans_Max_DD_Spring <- moran.test(fin.eu3$Max_DD_Spring, listw = w) %>% print()
morans_Mean_T_Summer <- moran.test(fin.eu3$Mean_T_Summer, listw = w) %>% print()
morans_Max_WD_Autumn <- moran.test(fin.eu3$Max_WD_Autumn, listw = w) %>% print()
morans_Max_DD_Autumn <- moran.test(fin.eu3$Max_DD_Autumn, listw = w) %>% print()
morans_gdp <- moran.test(fin.eu3$gdp, listw = w) %>% print()
morans_popdens <- moran.test(fin.eu3$pop_dens, listw = w) %>% print()

### Creeting a SpatialPointsDataFrame object (This is necessary because is the required format for the GWmodel library)
#spdf <- SpatialPointsDataFrame(coords, scaled_data)
#spdf_metric <- spTransform(spdf, CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")) # Transform to a metric CRS (e.g., UTM zone 33N)

#### Start from the last metric data frame - remove geometry column
tmp2 <- st_drop_geometry(scaled_data_meter) 

#### Set the CRS from the sf object
crs_proj4 <- st_crs(scaled_data_meter)$proj4string # Set the CRS from the sf object

#### Create a SPDF object
spdf <- SpatialPointsDataFrame(coords, tmp2, proj4string = CRS(crs_proj4))
str(spdf) # Check

#### Transform CRS to e.g., UTM (metric system)
spdf_metric <- spTransform(spdf, CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")) # transform the CRS to a metric system (e.g., UTM)
proj4string(spdf_metric) # Check

## Adjust Y using the spatial weights matrix
spdf_metric$Y_weighted <- W_matrix %*% spdf_metric$incidence
summary(spdf_metric$Y_weighted)

#spdf$Y_weighted <- W_matrix %*% spdf$incidence

# REGRESSION ANALYSIS ##########################################################
# GWR
## formula
 formula <- Y_weighted ~ Forest + Shrub + Urban + Crop + Water + Other + 
  Mean_P_Winter + Max_WD_Winter + Max_DD_Winter + Mean_T_Spring + Max_DD_Spring + 
  #Mean_T_Summer + 
  Max_WD_Autumn + Max_DD_Autumn + gdp + pop_dens

## Selecting the bandwidth
bw <- GWmodel::bw.gwr(formula = formula, 
                      data = spdf_metric, 
                      approach = "AIC", 
                      kernel = "bisquare",
                      adaptive = T)

## Checking for multicollinearity
### Local multicollinearity
lmc <- GWmodel::gwr.collin.diagno(formula = formula, 
                                 data = spdf_metric, 
                                 bw = bw,
                                 kernel = "bisquare",
                                 adaptive = T)
lmc

### Global multicollinearity
library(car)
lm <- lm(formula, spdf_metric)
gmc <- vif(lm)
gmc

# Fit GWR model
gwr_model <- gwr.basic(formula = formula,
                       data = spdf_metric, 
                       bw = bw,
                       kernel = "bisquare",
                       adaptive = T,
                       cv = T)

gwr_model
results <- gwr_model$SDF
min(results$Local_R2)
mean(results$Local_R2)
max(results$Local_R2)

# Quick plot R2
fin.eu3$r2 <- results$Local_R2 
fig.3.1 <- ggplot(fin.eu3) +
  geom_sf(aes(fill = r2, geometry = geometry))+
  scale_fill_viridis_c(na.value = "white", limits = c(0.75, 1), oob = scales::squish, direction = -1) +
  theme_ipsum() +
  labs(title = "Model diagnostics",
       subtitle="A) Coefficient of determination (R2), by NUTS3",
       fill="Local R2") +
  theme(plot.margin=unit(c(1,-2,-0.5,0),"cm"),
        legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=14),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=18))

fig.3.1

#
#
#
# Descriptive results (for the manuscript) #####################################
coef <- as_tibble(gwr_model$SDF) %>% 
  dplyr::select(#Forest, 
    Shrub, Urban, #Crop, 
    Water, #Other, 
                Mean_P_Winter, #Max_WD_Winter, 
    Max_DD_Winter, 
                Mean_T_Spring, Max_DD_Spring, 
                #Mean_T_Summer, 
                Max_WD_Autumn, Max_DD_Autumn, 
                gdp, pop_dens) %>%
  mutate(across(.cols = everything(), .fns= abs, .names = "abs_{col}"),
         geometry = scaled_data_meter$geometry)

coef.plot <- coef %>%
  mutate(Max_Coef = pmax(#abs_Forest, 
    abs_Shrub, abs_Urban, #abs_Crop, 
    abs_Water, #abs_Other, 
                         abs_Mean_P_Winter, #abs_Max_WD_Winter, 
    abs_Max_DD_Winter, 
                         abs_Mean_T_Spring, abs_Max_DD_Spring, 
                         #abs_Mean_T_Summer, 
                         abs_Max_WD_Autumn, abs_Max_DD_Autumn, 
                         abs_gdp, abs_pop_dens),
         Max_Origin = case_when(
           #Max_Coef == abs_Forest ~ "Forest",
           Max_Coef == abs_Shrub ~ "Shrub",
           Max_Coef == abs_Urban ~ "Urban",
           #Max_Coef == abs_Crop ~ "Crop",
           Max_Coef == abs_Water ~ "Water",
           #Max_Coef == abs_Other ~ "Other",
           Max_Coef == abs_gdp ~ "GDP",
           Max_Coef == abs_pop_dens ~ "Pop_dens",
           Max_Coef == abs_Mean_P_Winter ~ "Mean_P_Winter",
           #Max_Coef == abs_Max_WD_Winter ~ "Max_WD_Winter",
           Max_Coef == abs_Max_DD_Winter ~ "Max_DD_Winter",
           Max_Coef == abs_Mean_T_Spring ~ "Mean_T_Spring",
           Max_Coef == abs_Max_DD_Spring ~ "Max_DD_Spring",
           #Max_Coef == abs_Mean_T_Summer ~ "Mean_T_Summer",
           Max_Coef == abs_Max_WD_Autumn ~ "Max_WD_Autumn",
           Max_Coef == abs_Max_DD_Autumn ~ "Max_DD_Autumn")) %>%
  as_tibble()

#
#
#
# Figures (for the manuscript) #####################################
#
# Figure 1 - Outcome of interest
## Fig.1.1 - Map of incidence
fig1.tmp <- readRDS("20240627_APES FINAL_metric.rds") |>
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") |>
  dplyr::select(incidence, population, geometry) |>
  group_by(geometry) |>
  mutate(incidence_pck = (incidence * 100000),
         incidence_cat = case_when(is.na(incidence_pck) ~ NA,
                                   incidence_pck >= 100 ~ "A",
                                   incidence_pck < 100 & incidence_pck >= 10 ~ "B",
                                   incidence_pck < 10 & incidence_pck >= 1 ~ "C",
                                   incidence_pck < 1 & incidence_pck >= 0.1 ~ "D",
                                   incidence_pck < 0.1 & incidence_pck >= 0.01 ~ "E",
                                   incidence_pck < 0.01 ~ "FF"),
                  incidence_pct = (incidence * 100 / population) * 100000) |>
  as.data.frame() 
fig1.tmp$incidence_cat <- as.factor(fig1.tmp$incidence_cat)
table(fig1.tmp$incidence_cat)

#color.incidence <- c(FF = "#FA8072", E = "#CD5C5C", D ="#ED2939", C = "#BF0A30", B = "#B80F0A", A = "#7c0A02")   #Scale of color
#color.incidence <- c(FF = "#ffffc2", E = "#d16002",  D = "#d16002", C ="#c71585", B = "#702963", A = "#2c041d")  #Ema's palette
color.incidence <- c(FF = "#e7d5c5", E = "#ffbb55",  D = "#ffbb55", C ="#c44e4f", B = "#143c5d", A = "#011c39")  #Mindfull palettes

check <- fig1.tmp |> filter(incidence_cat != "FF")
fig.1 <- ggplot(fig1.tmp) +
  geom_sf(aes(fill = incidence_cat, geometry = geometry), linewidth = 0.05) + 
  geom_sf(data = check, aes(fill = incidence_cat, geometry = geometry), linewidth = 0.05) +
  theme_ipsum() + 
  labs(title = "Total incidence of West Nile virus human cases ",
       subtitle = "Expressed as total cases/10,000, at NUTS3 level in Continental Europe (2005-2019)",
       #subtitle = "A) Expressed as total cases/10,000, at NUTS3 level in Continental Europe (2005-2019)",
       fill = "Incidence") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=14),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=18)) +
  scale_fill_manual(values = color.incidence, labels = c(">=100", "100-10", "10-1", "1-0.1", #"0.1-0.01", # There are no values within the 0-1-0.01 category
                                                         "<0.1"), na.value = "white") 

fig.1

## Fig.1.2 - Incidence distribution
fig.1.2 <- ggplot(fig1.tmp) +
  geom_boxplot(aes(y = incidence_pck, fill = incidence_pck), linewidth = 0.15, alpha = 0.5) +
  scale_fill_gradient(low = "#FA8072", high = "#7c0A02", na.value = "white") +
  theme_ipsum() +
  xlim(-0.2, 0.2) +
  labs(title = " ",
       subtitle = "B) Distribution of the incidence/100,000 inhabitants") +
  theme(plot.subtitle = element_text(size=14)) 


fig.1.2

#fig.1 <- grid.arrange(fig.1, fig.1.2, ncol = 2, heights = c(2.5,1), widths = c(8,1))
ggsave(fig.1, filename = "Figure1_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")
ggsave(fig.1.2, filename = "Figure1.2_APES.png", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

fig.3bis <- grid.arrange(fig.3.1, fig.3.3, ncol = 2, heights = c(5,1), widths = c(2.5,0.85))
ggsave(fig.3bis, filename = "Figure3bis_APES.png", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

#
# Figure 2 - Covariates
## Fig.2.1 - 
avg_var <- c("Forest", "Shrub", "Urban", "Crop", "Water", "Other", 
             "Mean_P_Winter", "Max_WD_Winter", "Max_DD_Winter", 
             "Mean_T_Spring", "Max_DD_Spring", 
             "Mean_T_Summer",
             "Max_WD_Autumn", "Max_DD_Autumn", 
             "gdp", "pop_dens")

fig2.tmp <- readRDS("20240627_APES FINAL_metric.rds") |>
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") |>
  dplyr::select(incidence, Forest, Shrub, Urban, Crop, Water, Other, 
                Mean_T_Winter, Mean_T_Spring, Mean_T_Summer, Mean_T_Autumn,
                Mean_P_Winter, Mean_P_Spring, Mean_P_Summer, Mean_P_Autumn,
                Max_WD_Winter, Max_DD_Winter, Max_DD_Spring, Max_WD_Autumn, 
                Max_DD_Autumn, 
                gdp, pop_dens, geometry, NUTS_CODE, RecordType) |>
  group_by(geometry, NUTS_CODE, RecordType) |>
  dplyr::summarise(incidence = sum(incidence, na.rm = TRUE), 
                   across(all_of(avg_var), mean, na.rm = TRUE)) |>
  mutate(cent = sf::st_centroid(geometry),
         Max_Cov = pmax(Forest, Shrub, Urban, Crop, Water, Other),
         Max_Origin = case_when(
           Max_Cov == Forest ~ "Forest",
           Max_Cov == Shrub ~ "Shrub",
           Max_Cov == Urban ~ "Urban",
           Max_Cov == Crop ~ "Crop",
           Max_Cov == Water ~ "Water",
           Max_Cov == Other ~ "Other")) |> 
  na.omit() |>
  as_tibble()

group.colors <- c(Water = "#333BFF", Crop= "#E3DB71", Other ="#9633FF", 
                  Urban = "#ffa500", Forest = "#006400", Shrub = "#90EE90")

fig2.1 <- ggplot(fig2.tmp, aes(x = log(Forest))) +
  geom_density(fill = "#006400",  linewidth = 0.15, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Crop)), fill = "#E3db71", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Shrub)), fill = "#90EE90", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Urban)), fill = "#ffa500", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Water)), fill = "#333Bff", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Other)), fill = "#9633FF", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  theme_ipsum()+
  labs(title = "Land cover, demographic, and climatic covariates", 
       subtitle = "(A) Land cover distribution",
       x = "Log(Value)", y = "Density",color = "Land cover") +
  theme(plot.title = element_text(size=28),
        plot.subtitle = element_text(size=22))

fig2.1

fig2.2 <- ggplot(fig2.tmp) +
  geom_sf(aes(fill=Max_Origin, geometry=geometry), linewidth = 0.05, alpha = 0.75) + 
  theme_ipsum() +
  scale_fill_manual(values=group.colors) +
  labs(title = " ",
       subtitle="B) Main land cover category, by NUTS3",
       fill="Land cover") +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=18),
        plot.subtitle = element_text(size=22))

fig2.2

fig2.tmp_latlon <- fig2.tmp |>
  dplyr::mutate(lon = sf::st_coordinates(cent)[,1],
                lat = sf::st_coordinates(cent)[,2])

fig2.3 <- ggplot() +
  geom_sf(data = fig2.tmp_latlon, aes(geometry = geometry), fill = "white", linewidth = 0.05) + 
  geom_point(data = fig2.tmp_latlon, aes(x = lon, y = lat, size = log(gdp)), color = "#333Bff", 
             position = position_nudge(x = 0.25, y = 0), alpha = 0.75) +
  geom_point(data = fig2.tmp_latlon, aes(x = lon, y = lat, size = log(pop_dens)), color = "#ffa500",
             position = position_nudge(x = -0.25, y = 0), alpha = 0.75) +
  theme_ipsum() + 
  labs(title = "",
       subtitle = "C) Population density and GDP, by NUTS3",
       size = "Log of GDP (Yellow), Pop Dens (Blue)") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=18),
        plot.subtitle = element_text(size=22))
fig2.3

#fig2.t <- grid.arrange(grobs = list(fig2.1, fig2.2, fig2.3), ncol = 1, heights = c(1))
#ggsave(fig2, filename = "Figure2_APES.png", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

ggsave(fig2.1, filename = "Figure2.1_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")
ggsave(fig2.2, filename = "Figure2.2_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")
ggsave(fig2.3, filename = "Figure2.3_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

temp_long <- fig2.tmp %>%
  dplyr::select(Mean_T_Winter, Mean_T_Spring, Mean_T_Summer, Mean_T_Autumn) %>%
  rename(Winter = Mean_T_Winter, Spring = Mean_T_Spring, Summer = Mean_T_Summer, Autumn = Mean_T_Autumn) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value))

rain_long <- fin.eu.scale %>%
  select(Mean_P_Winter, Mean_P_Spring, Mean_P_Summer, Mean_P_Autumn) %>%
  rename(Winter = Mean_P_Winter, Spring = Mean_P_Spring, Summer = Mean_P_Summer, Autumn = Mean_P_Autumn) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value))
## Fig.2.2 - 

## Fig.2.3 - 

#
# Figure 3
## Fig.3.2 - Residuals distribution
# Fig.3.3 - Residuals distribution
resid <- gwr_model$lm$residuals %>% as_tibble()
fig.3.3 <- ggplot(resid, aes(y = value)) +
  geom_density(fill = "#1E90FF",  linewidth = 0.05, alpha = 0.5) +
  theme_ipsum() +
  #theme(plot.margin=unit(c(1,2,-0.5,-4.5),"cm")) + 
  labs(title = " ",
       subtitle = "B) Residuals' distribution") +
  theme(plot.subtitle = element_text(size=18))

fig.3.3
fig3 <- grid.arrange(fig.3.1, fig.3.3, ncol = 2, heights = c(5,1), widths = c(2.5,0.85))
ggsave(fig3, filename = "Figure3_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

#
# Fig.4 - Local Coefficients
coef <- as_tibble(gwr_model$SDF) %>% 
  dplyr::select(Shrub, Urban, Water, Mean_P_Winter, Max_DD_Winter, Mean_T_Spring, Max_DD_Spring,
    Max_WD_Autumn, Max_DD_Autumn, gdp, pop_dens) %>%
  mutate(across(.cols = everything(), .fns= abs, .names = "abs_{col}"),
         geometry = scaled_data_meter$geometry)

coef.plot <- coef %>%
  mutate(Max_Coef = pmax(abs_Shrub, abs_Urban, abs_Water,  
    abs_Mean_P_Winter, abs_Max_DD_Winter, abs_Mean_T_Spring, abs_Max_DD_Spring, abs_Max_WD_Autumn, abs_Max_DD_Autumn, 
    abs_gdp, abs_pop_dens),
    Max_Origin = case_when(
      Max_Coef == abs_Shrub ~ "Shrub",
      Max_Coef == abs_Urban ~ "Urban",
      Max_Coef == abs_Water ~ "Water",
      Max_Coef == abs_gdp ~ "gdp",
      Max_Coef == abs_pop_dens ~ "pop_dens",
      Max_Coef == abs_Mean_P_Winter ~ "Mean_P_Winter",
      Max_Coef == abs_Max_DD_Winter ~ "Max_DD_Winter",
      Max_Coef == abs_Mean_T_Spring ~ "Mean_T_Spring",
      Max_Coef == abs_Max_DD_Spring ~ "Max_DD_Spring",
      Max_Coef == abs_Max_WD_Autumn ~ "Max_WD_Autumn",
      Max_Coef == abs_Max_DD_Autumn ~ "Max_DD_Autumn")) %>%
  as_tibble()

## Fig. 4.1 - Map with the highest coefficient
group.colors <- c(Shrub = "#2E8B57", Urban = "#355E3B", Water ="#008080",  
                  gdp = "#ffa700", pop_dens = "#ffD300",
                  Mean_P_Winter = "#CAF0F8", Max_WD_Winter = "#90E0EF", Max_DD_Winter = "#7adaec",
                  Mean_T_Spring ="#48CAE4", Max_DD_Spring = "#00B4D8", Max_WD_Autumn = "#0077B6", Max_DD_Autumn = "#023E8A")

addmargins(table(coef.plot$Max_Origin))
coef.plot$Max_Origin <- coef.plot$Max_Origin |> 
  factor(levels = c("Shrub", "Urban", "Water", "gdp", "pop_dens", "Mean_P_Winter", "Max_WD_Winter", "Mean_T_Spring", "Max_DD_Spring", "Max_WD_Autumn", "Max_DD_Autumn"))

fig.4.1 <- ggplot() +
  geom_sf(data = coef.plot, aes(fill = Max_Origin, geometry = geometry, alpha = log(Max_Coef)), linewidth = 0.05) +
  theme_ipsum() +
  scale_fill_manual(values = group.colors, na.value = "white") +
  labs(title = "Coefficient estimates",
       subtitle = "A) Main coefficient estimate, by NUTS3",
       fill = "Main coefficient",
       alpha = "Coefficient value") +
  guides(fill = guide_legend(nrow = 2),
         alpha = guide_legend(nrow = 1),)+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=14),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=18))


fig.4.1
ggsave(fig.4.1, filename = "Figure4.1_APES.png", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")
#
#
#
## Fig. 4.2 - Boxplot
coef_long <- coef %>%
  dplyr::select(Shrub, Urban, Water, 
                Mean_P_Winter,Max_DD_Winter, 
                Mean_T_Spring, Max_DD_Spring, 
                Max_WD_Autumn, Max_DD_Autumn, 
                gdp, pop_dens) %>%
  pivot_longer(cols = everything(),
               names_to = "Variable", 
               values_to = "Max_Coef")

 coef_long$Variable <- coef_long$Variable |> 
  factor(levels = c("Shrub", "Urban", "Water", "gdp", "pop_dens", "Mean_P_Winter", "Max_DD_Winter", "Mean_T_Spring", "Max_DD_Spring", "Max_WD_Autumn", "Max_DD_Autumn"))
# Create the plot
fig.4.2 <- ggplot(coef_long, aes(x = Variable, y = Max_Coef)) +
  geom_boxplot(aes(fill = Variable), linewidth = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = group.colors, na.value = "white") +
  theme_ipsum() +
  scale_y_continuous(breaks = seq(-10, 10, by = 2),limits = c(-10, 10)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(title = " ",       
       subtitle = "B) Coefficient distribution",
       x = " ",
       y = " ") +
  theme(legend.position = "none",
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size=14),
        plot.subtitle = element_text(size=18))
  
   #+
#theme(plot.margin=unit(c(1,2,-0.5,-2.5),"cm"))

fig.4.2
fig.4 <- grid.arrange(fig.4.1, fig.4.2, ncol = 2, heights = c(5,1), widths = c(2.2,1))
ggsave(fig.4, filename = "Figure4_APES.tiff", bg = "white", dpi=300, width = 22, height = 12.5, units = "in")

#
#
#
# WORK IN PROGRESS ###################################################
#
# Fig.2 ##############################
fig2.tmp <- readRDS("20240627_APES FINAL_metric.rds") |>
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") |>
  dplyr::select(incidence, Forest, Shrub, Urban, Crop, Water, Other, 
                Mean_T_Winter, Mean_T_Spring, Mean_T_Summer, Mean_T_Autumn,
                Mean_P_Winter, Mean_P_Spring, Mean_P_Summer, Mean_P_Autumn,
                Max_WD_Winter, Max_DD_Winter, Max_DD_Spring, Max_WD_Autumn, 
                Max_DD_Autumn, 
                gdp, pop_dens, geometry, NUTS_CODE, RecordType) |>
  group_by(geometry, NUTS_CODE, RecordType) |>
  dplyr::summarise(incidence = sum(incidence, na.rm = TRUE), 
                   across(all_of(avg_var), mean, na.rm = TRUE)) |>
  mutate(cent = sf::st_centroid(geometry),
         Max_Cov = pmax(Forest, Shrub, Urban, Crop, Water, Other),
         Max_Origin = case_when(
           Max_Cov == Forest ~ "Forest",
           Max_Cov == Shrub ~ "Shrub",
           Max_Cov == Urban ~ "Urban",
           Max_Cov == Crop ~ "Crop",
           Max_Cov == Water ~ "Water",
           Max_Cov == Other ~ "Other")) |> 
  na.omit() |>
  as_tibble()

group.colors <- c(Water = "#333BFF", Crop= "#E3DB71", Other ="#9633FF", 
                  Urban = "#ffa500", Forest = "#006400", Shrub = "#90EE90")

fig.2.1 <- ggplot(fig2.tmp, aes(x = log(Forest))) +
  geom_density(fill = "#006400",  linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Crop)), fill = "#E3db71", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Shrub)), fill = "#90EE90", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Urban)), fill = "#ffa500", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Water)), fill = "#333Bff", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  geom_density(aes(x = log(Other)), fill = "#9633FF", linewidth = 0.05, alpha = 0.5) +  # Add density curve
  theme_ipsum()+
  labs(title = "Land cover, demographic, and climatic covariates", 
       subtitle = "(A) Land cover distribution",
       x = "Log(Value)", y = "Density",color = "Land cover")

fig.2.1

fig2.2 <- ggplot(fig2.tmp) +
  geom_sf(aes(fill=Max_Origin, geometry=geometry), linewidth = 0.05, alpha = 0.75) + 
  theme_ipsum() +
  scale_fill_manual(values=group.colors) +
  labs(title = " ",
       subtitle="B) Main land cover category, by NUTS3",
       fill="Land cover") 

fig2.2

fig2.tmp_latlon <- fig2.tmp |>
  dplyr::mutate(lon = sf::st_coordinates(cent)[,1],
                lat = sf::st_coordinates(cent)[,2])

fig2.3 <- ggplot() +
  geom_sf(data = fig2.tmp_latlon, aes(geometry = geometry), fill = "white", linewidth = 0.05) + 
  geom_point(data = fig2.tmp_latlon, aes(x = lon, y = lat, size = log(gdp)), color = "#333Bff", 
             position = position_nudge(x = 0.25, y = 0), alpha = 0.75) +
  geom_point(data = fig2.tmp_latlon, aes(x = lon, y = lat, size = log(pop_dens)), color = "#ffa500",
             position = position_nudge(x = -0.25, y = 0), alpha = 0.75) +
  theme_ipsum() + 
  labs(title = "",
       subtitle = "C) Population density and GDP, by NUTS3",
       alpha = "Pop dens",
       size = "GDP")

fig2.3

temp_long <- fig2.tmp %>%
  dplyr::select(Mean_T_Winter, Mean_T_Spring, Mean_T_Summer, Mean_T_Autumn) %>%
  rename(Winter = Mean_T_Winter, Spring = Mean_T_Spring, Summer = Mean_T_Summer, Autumn = Mean_T_Autumn) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value))

rain_long <- fin.eu.scale %>%
  select(Mean_P_Winter, Mean_P_Spring, Mean_P_Summer, Mean_P_Autumn) %>%
  rename(Winter = Mean_P_Winter, Spring = Mean_P_Spring, Summer = Mean_P_Summer, Autumn = Mean_P_Autumn) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Value))

#
# Fig.1 ##############################
fig1.tmp <- readRDS("20240627_APES FINAL_metric.rds") |>
  filter(!st_is_empty(geometry),
         RecordType != "ZIKV") |>
  dplyr::select(incidence, population, geometry) |>
  group_by(geometry) |>
  mutate(incidence_pck = (incidence * 100000),
         incidence_cat = case_when(is.na(incidence_pck) ~ NA,
                                   incidence_pck >= 100 ~ "A",
                                   incidence_pck < 100 & incidence_pck >= 10 ~ "B",
                                   incidence_pck < 10 & incidence_pck >= 1 ~ "C",
                                   incidence_pck < 1 & incidence_pck >= 0.1 ~ "D",
                                   incidence_pck < 0.1 & incidence_pck >= 0.01 ~ "E",
                                   incidence_pck < 0.01 ~ "F",
                                   TRUE ~ NA_character_),
         incidence_test = ifelse(incidence_pck >= 100, "A", "B"),
        incidence_pct = (incidence * 100 / population) * 100000) |>
  as.data.frame() |>
  st_as_sf()

table(fig1.tmp$incidence_cat)
fig1.tmp$incidence_cat <- as.factor(fig1.tmp$incidence_cat)  
table(fig1.tmp$incidence_test)
fig1.tmp$incidence_test <- as.factor(fig1.tmp$incidence_test)  
check <- fig1.tmp |>
  filter(incidence_test == "A")
  
st_geometry(fig1.tmp$geometry)

fig.1 <- ggplot(fig1.tmp) +
  geom_sf(aes(fill = incidence_pck, geometry = geometry), linewidth = 0.05))
          
          , geometry = geometry))



  theme_ipsum() +
  scale_fill_viridis_d() 

+
  labs(title = "Coefficient estimates", 
       subtitle = "A) Main Coefficient Estimate, by NUTS3", 
       fill = "Coef") +
  theme(legend.position = "none") +
  guides(alpha = "none") + # Remove the alpha legend
  theme(plot.margin=unit(c(1,0.5,-0.5,1),"cm")) 

fig1.2 <- ggplot(fig1.tmp, aes(y = incidence_pck)) +
  geom_boxplot(fill = "#1E90FF",  linewidth = 0.05, alpha = 0.5) +
  scale_y_continuous(limits=c(0.001,100)) +
  theme_ipsum() 


fig1.2
+
  #theme(plot.margin=unit(c(1,2,-0.5,-4.5),"cm")) + 
  labs(title = " ",
       subtitle = "B) Residuals' distribution") 
       
       ,
         subtitle = "B) Coefficient distribution",
         x = " ",
         y = " ") +
    theme(legend.position = "right") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_y_continuous(breaks = seq(-50, 50, by = 25),limits = c(-50, 50)) #+
  #theme(plot.margin=unit(c(1,2,-0.5,-2.5),"cm"))
  
  
library(treemap)
data <- data.frame(
  category = c(rep("Land cover", 5), rep("Clima", 6), rep("Socio-demographic", 2)),
  subcategory = c("Shrub", "Urban", "Crop", "Water",
               "Other", "Mean_P_Winter", "Max_DD_Winter", "Mean_T_Spring",
               "Max_DD_Spring", "Mean_T_Summer", "Max_DD_Autumn", "GDP", "Pop_dens"),
  value = c(195, 3, 6, 3, 20, 24, 4, 181, 12, 572, 8, 78, 28)
)

ggplot(data, aes(fill = subcategory, x = value)) +
  geom_bar(position = "dodge")

ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity")
treemap(data,
        index = c("category", "subcategory"), vSize = "value", type = "categorical",
        fontsize.labels = c(15, 12), # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels = c("white", "orange"), # Color of labels
        fontface.labels = c(2, 1), # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels = c("transparent"), # Background color of labels
        align.labels = list(
          c("center", "center"),
          c("right", "bottom")
        ), # Where to place labels in the rectangle?
        overlap.labels = 0.5, # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels = F, # If true, labels are bigger when rectangle is bigger.
)
png(filename="tree.png",width=1150, height=800)
treemap(data,
        index = c("category", "subcategory"),  # Categories and subcategories
        vSize = "value",                       # Size of the tiles
        vColor = "category",                   # Color by category
        type = "categorical",                  # Type of color scale
        palette = "Set2",
        #title = "Most relevant factor, by NUTS",
        fontsize.labels = c(18, 12),             # Font size for labels
        fontcolor.labels = "black",              # Font color for labels
        fontface.labels = c(2, 1),             # Font face for labels (bold for categories)
        bg.labels = c("transparent"), # Background color of labels
        align.labels = list(c("center", "center"), c("right", "bottom")), # Label alignment
        border.col = c("black", "white"),                  # Border color for the tiles
        border.lwds = c(0.3, 0.2),                # Border width for the tiles
        inflate.labels = F
)

install.packages("treemapify")
library(treemapify)

#
# Fig.1 ##############################