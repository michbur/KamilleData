#install.packages("ghit")
#ghit::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
library(dplyr)
library(tabulizer)
library(reshape2)
system("pdftk ft.pdf cat 7-12 output tab.pdf")
pdf_tabs <- extract_areas("tab.pdf")
save(pdf_tabs, file = "./results/pdf_tabs.RData")

part1 <- lapply(pdf_tabs[-length(pdf_tabs)], function(i) {
  genera_name <- strsplit(i[, 1], " ") %>% 
    lapply(first) %>% 
    unlist(use.names = FALSE)
  
  subst <- i[, 2]
  data.frame(genus = genera_name, subs = subst, stringsAsFactors = FALSE) %>% 
    filter(subst != "", genus != "Genus") 
}) %>% 
  do.call(rbind, .)

part2 <- local({
  genera_name <- pdf_tabs[[length(pdf_tabs)]][, 1] %>% 
    strsplit(" ") %>% 
    lapply(first) %>% 
    unlist(use.names = FALSE)
  
  subst <- pdf_tabs[[length(pdf_tabs)]][, 1] %>% 
    strsplit(") ") %>% 
    lapply(last) %>% 
    unlist(use.names = FALSE) 
  
  data.frame(genus = genera_name, subs = subst, stringsAsFactors = FALSE) %>% 
    filter(!(genus %in% c("Genus", "Order", "Class", "areferences"))) 
})

subs_dat <- rbind(part1, part2) %>% 
  filter(!grepl("[0-9]", x = genus)) %>% 
  mutate(subs = ifelse(subs == "Stenotrophomonas Palleroni and Bradbury Pyrene", "Pyrene", subs)) %>% 
  mutate(subs_rel = subs %in% c("Benzo[a]pyrene", "Crude oil", "Pyrene"))



otu_meta <- read.table("./data/otus_tax_assignments.txt", sep = "\t", header = TRUE)

otu_names <- select(otu_meta, OTUID, genus) %>% 
  filter(!is.na(genus)) %>% 
  select(genus) %>% 
  unlist(use.names = FALSE) %>% 
  as.character()

# check for misnamed genera
setdiff(filter(subs_dat, subs_rel)[["genus"]],
        otu_names) %>%
  lapply(function(i) {
    list(name = i,
         otu_names = unique(otu_names[grep(pattern = substr(i, 2, 6), x = otu_names)]))
    })


otu_full <- read.table("./data/otu_table_norm_samples.txt", sep = "\t", header = TRUE) %>% 
  melt %>% 
  inner_join(otu_meta %>% 
               select(OTUID, genus) %>% 
               left_join(filter(subs_dat, subs_rel) %>% select(-subs))
  ) %>% 
  mutate(ap = ifelse(is.na(genus), NA, ifelse(is.na(subs_rel), FALSE, TRUE))) %>% 
  select(-subs_rel)
  

otu_frac <- otu_full %>% 
  group_by(variable, ap) %>% 
  summarise(count = sum(value)) %>% 
  mutate(frac = count/sum(count)) %>% 
  ungroup() %>% 
  mutate(ap = factor(ap, exclude = NULL, labels = c("non-processing", "processing", "unknown"))) %>% 
  mutate(ap = factor(ap, levels = c("unknown", "non-processing", "processing")))

# library(ggplot2)
# ggplot(otu_frac, aes(x = variable, fill = ap, y = frac)) +
#   geom_bar(stat = "identity", position = "stack")

otu_frac %>% 
  filter(ap != "unknown") %>% 
  group_by(variable) %>% 
  mutate(frac = count/sum(count)) %>% 
  ungroup() %>% 
  filter(ap == "processing") %>% 
  mutate(SampleID = gsub(".", "-", variable, fixed = TRUE)) %>% 
  select(SampleID, frac) %>% 
  write.csv("./results/sample_processing.csv", row.names = FALSE)
 
write.csv(otu_frac, "./results/full_processing.csv", row.names = FALSE)
