library(dplyr)
library(ggplot2)

dat <- read.table("data/KamilleData", header = TRUE)

dat_g <- filter(dat, Origin == "Greenland") %>% 
  droplevels()

dat_s <- filter(dat, Origin == "Svalbard") %>% 
  droplevels()

qqnorm(dat_g[["Shannon"]])
qqline(dat_g[["Shannon"]])

summary(lm(Shannon ~ Species * Treatment * TemperatureC, data = dat_g))

summary(lm(Shannon ~ Treatment * SamplingDate, data = dat_s))
