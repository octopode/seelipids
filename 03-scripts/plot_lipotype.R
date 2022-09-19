library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))

# load the parsed data
file_tidy = here("02-tidydata", "lipotype_tidy.tsv")
# skip if parsing was just done and ltypdata_long is already in the global environment
ltypdata_long = read_tsv(file_tidy)

# this file annotates unique extract IDs (`eid`) with metadata (`sp`, `treatment`, etc)
file_meta = here("01-rawdata", "lipotype_metadata.tsv")
# load metadata
metadata = file_meta %>% read_tsv()

# join metadata to lipid data
ltypdata = ltypdata_long %>%
  ungroup() %>% 
  left_join(metadata, by = "eid") %>% 
  # then make compound ID a factor so it plots in a uniform order
  mutate(
    # put headgroups in the order specified in the color palette,
    # which can be changed in seelipids_helpers.R
    class = factor(class, levels = names(chroma_cl)),
    # then order things by total chain length and unsaturation
    id = factor(
      id,
      levels = cur_data() %>%
        distinct(id, class, carbon, dbonds) %>%
        arrange(class, carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>% 
  # sort the rows for easy viewing
  arrange(eid, id)

# store it
ltypdata %>% write_tsv(here("02-tidydata", "lipotype_wmeta.tsv"))

#### PLOTS
### L. kluyveri pressure incubation experiments
### 2 strains x 2 pressures x 3 bio. replicates

## PLOT1: all lipid species, plotted individually
ltypdata %>% 
  # these are the only L. kluyveri in the dataset
  filter(sp == "Lach_kluy") %>% 
  # order the strains
  mutate(strain = factor(strain, levels = c("FM628", "fad2"))) %>% 
  gg_headgp(
    darkmode = FALSE,
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
    title = "L. kluyveri total lipids",
    x = "Biological replicate",
    y = "Mole fraction total lipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "pressureyeast_tot_lipids.pdf"), width = 10, height = 8)
ggsave(here("05-png", "pressureyeast_tot_lipids.png"), width = 10, height = 8)

## PLOT2: phospholipids, plotted individually
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
    darkmode = FALSE,
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
    title = "L. kluyveri phospholipids",
    x = "Biological replicate",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "pressureyeast_phospholipids.pdf"), width = 10, height = 8)
ggsave(here("05-png", "pressureyeast_phospholipids.png"), width = 10, height = 8)

## PLOT3: phospholipids, lumped by class and plotted individually
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
  # sum compounds in each class
  group_by(press, strain, eid, rep, class) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = FALSE,
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
    title = "L. kluyveri phospholipids",
    x = "Biological replicate",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "pressureyeast_phospholipids_lumped.pdf"), width = 10, height = 8)
ggsave(here("05-png", "pressureyeast_phospholipids_lumped.png"), width = 10, height = 8)

## PLOT4: phospholipids only, lumped by class, averaged by strain x pressure,
## laid out a little differently
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
  # sum compounds in each class
  group_by(eid, class) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # average biological replicates
  group_by(strain, press, class) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = factor(press),
    y = frac_molar,
    fill = class
  ) +
  facet_grid(
    cols = vars(strain)
  ) +
  labs(
    title = "L. kluyveri phospholipids",
    x = "Pressure (bar)",
    y = "Mean mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "pressureyeast_phospholipids_strainmean.pdf"), width = 6, height = 4)
ggsave(here("05-png", "pressureyeast_phospholipids_strainmean.png"), width = 6, height = 4)

## PLOT5: PL acyl carbon distributions within a few headgroup classes
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
  # sum compounds by saturation level
  group_by(eid, dbonds) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # average biological replicates
  group_by(strain, press, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  group_by(strain, press) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 1.0,
    x = dbonds,
    y = frac_molar,
    fill = strain,
    cols = vars(strain),
    rows = vars(press)
  ) +
  scale_fill_manual(
    values = brewer.pal(3, "Set2")[1:2] %>% 
      setNames(c("FM628", "fad2"))
  ) +
  labs(
    title = "PL acyl double bond distribution",
    x = "Total acyl double bonds",
    y = "Mole fraction phospholipids"
  )
# save vector and raster images
ggsave(here("04-pdf", "pressureyeast_acyldbonds.pdf"), width = 6, height = 6)
ggsave(here("05-png", "pressureyeast_acyldbonds.png"), width = 6, height = 6)

### S. cerevisiae CL experiments
