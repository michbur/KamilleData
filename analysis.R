library(dplyr)
library(ggplot2)
library(broom)

dat <- read.table("data/KamilleData", header = TRUE)

# analysis Greenland ----------------------------------------------

dat_g <- filter(dat, Origin == "Greenland") %>% 
  droplevels() %>% 
  #group_by(SampleID) %>% 
  mutate(replic = unlist(lapply(strsplit(as.character(SampleID), "-"), last)),
         id = factor(as.numeric(factor(unlist(lapply(strsplit(as.character(SampleID), "-"), 
                                              function(i) paste0(i[-length(i)], collapse = ""))))))) %>% 
  group_by(id, Species, Treatment, TemperatureC) %>% 
  summarise(Shannon = median(Shannon)) %>% 
  ungroup %>% 
  mutate(#TemperatureC = TemperatureC > 2,
         Treatment = factor(Treatment, levels = c("Control", "Acetone_control", "1_nM_pyrene", 
                                                  "100_nM_pyrene", "100+_nM_pyrene"), ordered = TRUE))

#qqnorm(dat_s[["Shannon"]])
#qqline(dat_s[["Shannon"]])

m1 <- lm(Shannon ~ Species * TemperatureC * Treatment, data = dat_g)
m2 <- lm(Shannon ~ Species * TemperatureC + Treatment, data = dat_g)
anova(m1, m2)
summary(m1)
summary(m2)

# https://en.wikipedia.org/wiki/Akaike_information_criterion
# https://en.wikipedia.org/wiki/Bayesian_information_criterion
glance(m1)
glance(m2)

tidy(m2)

ggplot(dat_g, aes(x = Species, y = Shannon, color = Treatment, shape = as.factor(TemperatureC))) +
  geom_quasirandom(size = 5)

# analysis Svalbard ----------------------------------------------

#contrasts(dat_g[["Treatment"]]) <- contr.helmert 

dat_s <- filter(dat, Origin == "Svalbard") %>% 
  droplevels() %>% 
  mutate(replic = unlist(lapply(strsplit(as.character(SampleID), "-"), last)),
         id = factor(as.numeric(factor(unlist(lapply(strsplit(as.character(SampleID), "-"), 
                                                     function(i) paste0(i[-length(i)], collapse = ""))))))) %>% 
  #group_by(id, Species, Treatment, SamplingDate) %>% 
  #summarise(Shannon = median(Shannon)) %>% 
  ungroup() %>% 
  filter(Treatment != "Control_start") %>% 
  mutate(SamplingDate = factor(SamplingDate, labels = c("early spring", "late spring"))) %>% 
  droplevels() %>% 
  mutate(Treatment = relevel(Treatment, "Control_end"))

ggplot(dat_s, aes(x = SamplingDate, y = Shannon, color = Treatment)) +
  geom_quasirandom(size = 5) +
  theme_bw() 

m_s <- lm(Shannon ~ SamplingDate*Treatment, data = dat_s)

summary(m_s)
tidy(m_s)
glance(m_s)

interaction.plot(dat_s[["SamplingDate"]], dat_s[["Treatment"]], response = dat_s[["Shannon"]])
