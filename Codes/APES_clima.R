# Special characters {} [] ~

library(dplyr)
library(lubridate)
#clima05 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_2005.csv") %>% 
#  mutate(Years = as.Date(time),
#         temp = t2m - 273.15) %>%
#  select(-c(time, t2m)) %>%
#  as_tibble()

clima06 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2006.csv") %>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima06$Years) <- 2006

clima07 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2007.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima07$Years) <- 2007

clima08 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2008.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima08$Years) <- 2008

clima09 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2009.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima09$Years) <- 2009

clima10 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2010.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima10$Years) <- 2010

clima11 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2011.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima11$Years) <- 2011

clima12 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2012.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima12$Years) <- 2012

clima13 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2013.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima13$Years) <- 2013

clima14 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2014.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima14$Years) <- 2014

clima15 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2015.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima15$Years) <- 2015

clima16 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2016.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima16$Years) <- 2016

clima17 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2017.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima17$Years) <- 2017

clima18 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2018.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima18$Years) <- 2018

clima19 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2019.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima19$Years) <- 2019

clima20 <- read.csv("/eos/jeodpp/data/projects/APES/NUTS3_Climate/t2m_mean_2020.csv")%>%
  mutate(Years = as.Date(time),
         temp = t2m - 273.15) %>%
  select(-c(time, t2m)) %>%
  as_tibble()
year(clima20$Years) <- 2020

clima <- rbind(clima05, clima06, clima07, clima08, clima09, clima10, clima11, 
               clima12, clima13, clima14, clima15, clima16, clima17, clima18, 
               clima19, clima20) %>% 
    mutate(Year = lubridate::year(Years),
         Month = lubridate::month(Years),
         Day = lubridate::day(Years),
         Growth.day = ifelse(temp > 14.3, 1, 0)) %>% 
  rename(Date = Years) %>%
  group_by(NUTS_CODE, Year, Month)%>%
  mutate(Cons.GD = ave(Growth.day, cumsum(Growth.day == 0), FUN=cumsum)) %>%
  ungroup() %>%
  as_tibble()

 clima.month <- clima %>%
  group_by(NUTS_CODE, Year, Month) %>%
  summarise(mean.temp = mean(temp, na.rm = T),
            Growth.days.tot = sum(Growth.day, na.rm = T),
            Growth.days.max = max(Cons.GD)) %>% 
  as_tibble()

 clima.season <- clima.month %>%
   mutate(Season = case_when(Month == 1 | Month == 2 | Month == 3 ~ "Winter",
                           Month == 4 | Month == 5 | Month == 6 ~ "Spring",
                           Month == 7 | Month == 8 | Month == 9 ~ "Summer",
                           Month == 10 | Month == 11 | Month == 12 ~ "Autumn")) %>%
   group_by(NUTS_CODE, Year, Season) %>%
   summarise(mean.temp = mean(mean.temp, na.rm = T),
             Growth.days.tot = sum(Growth.days.tot, na.rm = T),
             Growth.days.max = max(Growth.days.max)) %>%
   as_tibble()
 
saveRDS(clima, file = "Clima.rds")
saveRDS(clima.month, file = "Clima_month.rds")
saveRDS(clima.season, file = "Clima_season.rds")
#
#
# reshape seasonal dataset
clima <- readRDS("Clima_season.rds")
win <- clima %>% 
  filter(Season == "Winter") %>% 
  rename(Mean_T_Winter = mean.temp,
         Tot_GD_Winter = Growth.days.tot,
         Max_GD_Winter = Growth.days.max,
         #Mean_Pre_Winter = Mean_Pre
  ) %>%
  select(-c(Season)) %>%
  as_tibble()

spr <- clima %>% 
  filter(Season == "Spring") %>% 
  rename(Mean_T_Spring = mean.temp,
         Tot_GD_Spring = Growth.days.tot,
         Max_GD_Spring = Growth.days.max,
         #Mean_Pre_Winter = Mean_Pre
  ) %>%
  select(-c(Season)) %>%
  as_tibble()

sum <- clima %>% 
  filter(Season == "Summer") %>% 
  rename(Mean_T_Summer = mean.temp,
         Tot_GD_Summer = Growth.days.tot,
         Max_GD_Summer = Growth.days.max,
         #Mean_Pre_Winter = Mean_Pre
  ) %>%
  select(-c(Season)) %>%
  as_tibble()

aut <- clima %>% 
  filter(Season == "Autumn") %>% 
  rename(Mean_T_Autumn = mean.temp,
         Tot_GD_Autumn = Growth.days.tot,
         Max_GD_Autumn = Growth.days.max,
         #Mean_Pre_Winter = Mean_Pre
  ) %>%
  select(-c(Season)) %>%
  as_tibble()


clima2 <- left_join(win, spr, by = c("Year", "NUTS_CODE"))
clima2 <- left_join(clima2, sum, by = c("Year", "NUTS_CODE"))
clima2 <- left_join(clima2, aut, by = c("Year", "NUTS_CODE"))
saveRDS(clima2, file = "clima2.rds")
