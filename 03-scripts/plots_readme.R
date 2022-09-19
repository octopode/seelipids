library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))
source(here("03-scripts", "parse_lipidmaps.R"))
source(here("03-scripts", "parse_lipotype.R"))

# a nice darkmode plot of PLs from LipoType
ltypdata %>% 
  # these are the only L. kluyveri in the dataset
  filter(sp == "Lach_kluy") %>% 
  # order the strains
  mutate(strain = factor(strain, levels = c("FM628", "fad2"))) %>% 
  # filter down to phospholipids
  filter(str_detect(class, 'P')) %>% 
  # and renormalize
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = TRUE,
    label_frac = 0.05, # only label compounds that are ≥5%
    # aesthetic mappings are passed straight thru
    x = rep,
    y = frac_molar,
    fill = class
  ) +
  facet_grid(
    rows = vars(strain),
    cols = vars(press)
  ) +
  labs(
    #title = "L. kluyveri phospholipids",
    x = "Biological replicate",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("05-png", "pressureyeast_phospholipids_dark.png"), width = 5, height = 4)

# a nice darkmode plot of chain length distros from 2 spp.
lmapdata %>% 
  # filter E. coli out of the dataset
  filter(!str_detect(sp, "Ecol")) %>% 
  # another filter for example
  filter(tissue %in% c("whole", "body")) %>% 
  # make sp and eid into ordered factors
  mutate(
    sp  = factor(sp,  levels = sp_order ),
    eid = factor(eid, levels = eid_order)
  ) %>% 
  # pick just a few headgroup classes
  filter(class %in% c("PC", "PE", "PS")) %>% 
  # and just 2 spp.
  filter(sp %in% c("Bath_fost", "Llyr_bent")) %>% 
  # average by species. 
  group_by(sp, class, id, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # add together all the compounds in each class with the same # acyl carbons
  # Comment this out to break each bar into sections for each individual compound.
  group_by(sp, class, carbon) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # normalize within each headgroup
  group_by(sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = TRUE,
    meanline = TRUE,
    # this shades the odd columns for contrast
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = class
  ) +
  labs(
    #title = "Ctenophore acyl carbon distribution by headgroup class and species",
    x = "Total acyl carbons",
    y = "Mole fraction of class"
  )
# save vector and raster images
ggsave(here("05-png", "somectenos_acylcarbons_dark.png"), width = 6, height = 4)
