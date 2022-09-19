library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))

# load the parsed data
file_tidy = here("02-tidydata", "lipidmaps_tidy.tsv")
# skip if parsing was just done and lmapdata_long is already in the global environment
lmapdata_long = read_tsv(file_tidy)

# this file annotates unique extract IDs (`eid`) with metadata (`sp`, `treatment`, etc)
file_meta = here("01-rawdata", "lipidmaps_metadata.tsv")
# load metadata
metadata = file_meta %>% read_tsv()

# join metadata to lipid data
lmapdata = lmapdata_long %>%
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
lmapdata %>% write_tsv(here("02-tidydata", "lipidmaps_wmeta.tsv"))

#### PLOTS

### ctenophores: wild animals annotated with continuous environmental metadata

## here's a little block of code that sorts the samples based on depth and temperature
## intended as an example of how samples can be ordered systematically
samps_sorted = lmapdata %>% 
  group_by(sp) %>% 
  # filter E. coli out of the dataset
  filter(!str_detect(sp, "Ecol")) %>% 
  mutate(
    # calculate species mean depth and temp
    across(c(depth, temp), mean, .names="{.col}_mean"),
    # assign species to temperature and depth slices
    cold = (temp_mean <= 10),
    deep = (depth_mean >= 200)
  ) %>% 
  ungroup() %>% 
  arrange(cold, -temp_mean, deep, depth_mean) %>% 
  # lock in species order as a factor
  mutate(sp = sp %>% factor(., unique(.))) %>% 
  # sort by individual temp and depth within sp
  arrange(sp, -temp, depth) %>% 
  group_by(sp, eid) %>% 
  summarize()

sp_order  = samps_sorted$sp %>% unique()
eid_order = samps_sorted$eid %>% unique()

#### if this is daunting, you can sort the samples yourself and hardcode the vector like so:
#sp_order = c("Boli_vitr", "Leuc_pulc", "Mert_angu", "Bath_fost", "Bero_cucu", "Cydi_blac", "Llyr_deep", "Llyr_bent", "Tjal_pink")

#eid_order = c("JWL0189", "JWL0190", "JWL0166", "JWL0167", "JWL0157", "JWL0159", 
#"JWL0171", "JWL0172", "JWL0173", "JWL0191", "JWL0174", "JWL0175", 
#"JWL0169", "JWL0170", "JWL0168", "JWL0176", "JWL0177", "JWL0178", 
#"JWL0179", "JWL0180", "JWL0181", "JWL0185", "JWL0188")

## PLOT1: all lipid species, plotted individually
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
  arrange(sp, eid) %>% 
  # make a new variable that shows sp next to eid
  unite("sp_eid", sp, eid, sep = ' ') %>% 
  mutate(sp_eid = sp_eid %>% factor(., levels = unique(.))) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "Ctenophore total lipids",
    x = "Species and sample ID",
    y = "Mole fraction total lipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "ctenos_tot_lipids.pdf"), width = 10, height = 8)
ggsave(here("05-png", "ctenos_tot_lipids.png"), width = 10, height = 8)

## PLOT2: phospholipids only, plotted individually
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
  arrange(sp, eid) %>% 
  # make a new variable that shows sp next to eid
  unite("sp_eid", sp, eid, sep = ' ') %>% 
  mutate(sp_eid = sp_eid %>% factor(., levels = unique(.))) %>% 
  # here's where we filter to phospholipids only
  group_by(sp_eid) %>%
  filter(str_detect(class, 'P')) %>% 
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "Ctenophore phospholipids",
    x = "Species and sample ID",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "ctenos_phospholipids.pdf"), width = 10, height = 8)
ggsave(here("05-png", "ctenos_phospholipids.png"), width = 10, height = 8)

## PLOT3: phospholipids only, lumped by class and averaged by species
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
  # filter to phospholipids only
  group_by(eid) %>%
  filter(str_detect(class, 'P')) %>% 
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  # now add together all the compounds in each class %>% 
  group_by(sp, eid, class) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # and average by species
  group_by(sp, class) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalize
  group_by(sp) %>% 
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = sp,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "Ctenophore phospholipids by species",
    x = "Species and sample ID",
    y = "Mean mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "ctenos_phospholipids_spmean.pdf"), width = 10, height = 8)
ggsave(here("05-png", "ctenos_phospholipids_spmean.png"), width = 10, height = 8)

## PLOT4: total acyl carbon distributions within a few headgroup classes
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
    darkmode = FALSE,
    meanline = TRUE,
    # this shades the odd columns for contrast
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "Ctenophore acyl carbon distribution by headgroup class and species",
    x = "Total acyl carbons",
    y = "Mole fraction of class"
  )
# save vector and raster images
ggsave(here("04-pdf", "ctenos_acylcarbons.pdf"), width = 10, height = 8)
ggsave(here("05-png", "ctenos_acylcarbons.png"), width = 10, height = 8)

## PLOT5: total double bond distributions within a few headgroup classes
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
  # average by species. 
  group_by(sp, class, id, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # add together all the compounds in each class with the same # acyl carbons
  # Comment this out to break each bar into sections for each individual compound.
  group_by(sp, class, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # normalize within each headgroup
  group_by(sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 1.0,
    x = dbonds,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "Ctenophore acyl double bond distribution by headgroup class and species",
    x = "Total acyl double bonds",
    y = "Mole fraction of class"
  )
# save vector and raster images
ggsave(here("04-pdf", "ctenos_acyldbonds.pdf"), width = 10, height = 8)
ggsave(here("05-png", "ctenos_acyldbonds.png"), width = 10, height = 8)

### E. coli: Balanced 2x2 experimental design - 2 strains x 2 culture pressures ='depths'

## PLOT6: phospholipids only, plotted individually in a 2x2 facet grid
lmapdata %>% 
  # filter NON-E. coli out of the dataset
  filter(str_detect(sp, "Ecol")) %>% 
  group_by(eid) %>% 
  arrange(sp, depth, eid) %>% 
  mutate(
    # make sp into an ordered factor so wild is on left and plpe on right
    sp  = factor(sp,  levels = c("Ecol_wild", "Ecol_plpe")),
    # and give each sample a simple replicate # for plotting purposes
    # note how they are grouped and sorted above
    # magic number 3 = number of bio replicates per treatment
    rep = cur_group_id() %% 3 + 1
  ) %>% 
  # here's where we filter to phospholipids only
  group_by(eid) %>%
  filter(str_detect(class, 'P')) %>% 
  # we will also remove P-PC, since the Lipid Maps folks determined
  # that it is a misannotation
  filter(class != "P-PC") %>% 
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # aesthetic mappings are passed straight thru
    x = rep,
    y = frac_molar,
    fill = class
  ) +
  # the facet grid to separate strains and pressures
  facet_grid(
    rows = vars(depth),
    cols = vars(sp)
  ) +
  labs(
    title = "E. coli phospholipids",
    x = " Biological replicate",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "Ecoli_phospholipids.pdf"), width = 6, height = 6)
ggsave(here("05-png", "Eccoli_phospholipids.png"), width = 6, height = 6)

## PLOT7: Aggregate acyl carbon count distributions by strain and pressure
lmapdata %>% 
  # filter NON-E. coli out of the dataset
  filter(str_detect(sp, "Ecol")) %>% 
  group_by(eid) %>% 
  arrange(sp, depth, eid) %>% 
  mutate(
    # make sp into an ordered factor so wild is on left and plpe on right
    sp  = factor(sp,  levels = c("Ecol_wild", "Ecol_plpe")),
    # and give each sample a simple replicate # for plotting purposes
    # note how they are grouped and sorted above
    # magic number 3 = number of bio replicates per treatment
    rep = cur_group_id() %% 3 + 1
  ) %>% 
  # here's where we filter to phospholipids only
  group_by(eid) %>%
  filter(str_detect(class, 'P')) %>% 
  # we will also remove P-PC, since the Lipid Maps folks determined
  # that it is a misannotation
  filter(class != "P-PC") %>% 
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  # average by species. 
  group_by(sp, depth, id, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # add together all the compounds in each class with the same # acyl carbons
  # Comment this out to break each bar into sections for each individual compound.
  group_by(sp, depth, carbon) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # normalize within each headgroup
  group_by(sp, depth) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    # this shades the odd columns for contrast
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = sp,
    # Here we assign custom colors to the strains
    # bc we want to override the default coding by headgroup class
    chroma = c(
      Ecol_wild = "lightgrey",
      Ecol_plpe = chroma_cl[["P-PE"]]
    ),
    # facet_wrap spec
    rows = vars(depth),
    cols = vars(sp)
  ) +
  labs(
    title = "E. coli phospholipid acyl carbon distribution by strain and pressure",
    x = "Total acyl carbons",
    y = "Mole fraction"
  )
# save vector and raster images
ggsave(here("04-pdf", "Ecoli_acylcarbons.pdf"), width = 10, height = 8)
ggsave(here("05-png", "Ecoli_acylcarbons.png"), width = 10, height = 8)

## PLOT8: Aggregate acyl double bond distributions by strain and pressure
lmapdata %>% 
  # filter NON-E. coli out of the dataset
  filter(str_detect(sp, "Ecol")) %>% 
  group_by(eid) %>% 
  arrange(sp, depth, eid) %>% 
  mutate(
    # make sp into an ordered factor so wild is on left and plpe on right
    sp  = factor(sp,  levels = c("Ecol_wild", "Ecol_plpe")),
    # and give each sample a simple replicate # for plotting purposes
    # note how they are grouped and sorted above
    # magic number 3 = number of bio replicates per treatment
    rep = cur_group_id() %% 3 + 1
  ) %>% 
  # here's where we filter to phospholipids only
  group_by(eid) %>%
  filter(str_detect(class, 'P')) %>% 
  # we will also remove PC and P-PC, since the Lipid Maps folks determined
  # that they are misannotations
  filter(!str_detect(class, "PC")) %>% 
  # the trace PUFAs are probably wrong too
  # comment to see them plotted
  filter(dbonds <= 2) %>% 
  # and renormalize to total phospholipids
  # this relies on the above group_by() call
  mutate(frac_molar = frac_molar / sum(frac_molar)) %>% 
  # average by species. 
  group_by(sp, depth, id, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # add together all the compounds in each class with the same # acyl carbons
  # Comment this out to break each bar into sections for each individual compound.
  group_by(sp, depth, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # normalize within each headgroup
  group_by(sp, depth) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    x = dbonds,
    y = frac_molar,
    fill = sp,
    # Here we assign custom colors to the strains
    # bc we want to override the default coding by headgroup class
    chroma = c(
      Ecol_wild = "lightgrey",
      Ecol_plpe = chroma_cl[["P-PE"]]
    ),
    # facet_wrap spec
    rows = vars(depth),
    cols = vars(sp)
  ) +
  labs(
    title = "E. coli phospholipid acyl double bond distribution by strain and pressure",
    x = "Total acyl double bonds",
    y = "Mole fraction"
  )
# save vector and raster images
ggsave(here("04-pdf", "Ecoli_acyldbonds.pdf"), width = 10, height = 8)
ggsave(here("05-png", "Ecoli_acyldbonds.png"), width = 10, height = 8)
