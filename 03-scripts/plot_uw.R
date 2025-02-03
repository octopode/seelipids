library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))

# load the parsed data
file_tidy = here("02-tidydata", "uw_tidy.tsv")
# skip if parsing was just done and uwdata_long is already in the global environment
uwdata_long = read_tsv(file_tidy)

# this file annotates unique extract IDs (`eid`) with metadata (`sp`, `treatment`, etc)
file_meta = here("01-rawdata", "uw_metadata.tsv")
# load metadata
metadata = file_meta %>% read_tsv()

# join metadata to lipid data
uwdata = uwdata_long %>%
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
  # sometimes it makes plotting easier by casting metadata
  # vars to factors with fixed order
  mutate(loc = loc %>% factor(., levels = c("rock", "intermed", "hardplace"))) %>% 
  # sort the rows for easy viewing
  arrange(eid, id)

# store it
uwdata %>% write_tsv(here("02-tidydata", "lipotype_wmeta.tsv"))

#### PLOTS
### I don't know what the samples represent, so I gave them fake species
### and locations as placeholders. Replace these with real metadata.

## PLOT1: all lipid species, plotted individually
uwdata %>% 
  mutate(
    # order the species in the desired plotting order
    sp = factor(sp, levels = c("unk2", "unk1")),
    # make EID a factor so it does not plot as a continuous variable
    eid = eid %>% as.factor()
  ) %>% 
  gg_headgp(
    darkmode = FALSE,
    # this can help declutter a complex plot by removing(!) low-abundance species
    # but if not left at 0, it should be set very low to avoid misrepresenting the data!
    #thres_draw = 0.002, 
    # how abundant a compound must be to get a text label
    label_frac = 0.01,
    # aesthetic mappings are passed straight thru
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  facet_grid(
    cols = vars(sp),
    #cols = vars(press)
    scales = "free_x" # so all factors do not plot in both panels
  ) +
  labs(
    title = "Headgroup composition by species and indl",
    x = "Extract ID",
    y = "Mole fraction total lipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_headgroupclasses.pdf"), width = 10, height = 8)
ggsave(here("05-png", "uw_headgroupclasses.png"), width = 10, height = 8)

## PLOT2: PL acyl carbon distributions
## (change filter to look at other classes)
uwdata %>% 
  # filter to select classes
  filter(class %in% c("CE")) %>% 
  # and renormalize
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # sum compounds by saturation level within each individual
  group_by(eid, class, dbonds) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and location
  group_by(sp, loc, class, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 1.0,
    x = dbonds,
    y = frac_molar,
    fill = sp,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(sp),
    rows = vars(loc)
  ) +
  scale_fill_manual(
    # specifies a color palette for sp
    values = brewer.pal(3, "Set2")[1:2] %>% 
      setNames(c("unk1", "unk2"))
  ) +
  labs(
    title = "CE acyl double bond distribution",
    x = "Total acyl double bonds",
    y = "Mole fraction phospholipids"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_acyldbonds.pdf"), width = 6, height = 6)
ggsave(here("05-png", "uw_acyldbonds.png"), width = 6, height = 6)

## PLOT3: Same as plot 2 but for chain lengths
uwdata %>% 
  # filter to select classes
  filter(class %in% c("CE")) %>% 
  # and renormalize
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # sum compounds by chain length within each individual
  group_by(eid, class, carbon) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and location
  group_by(sp, loc, class, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 0.4, # give lighter shading for odd-chain-length bars
    x = carbon,
    y = frac_molar,
    fill = sp,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(sp),
    rows = vars(loc)
  ) +
  scale_fill_manual(
    # specifies a color palette for sp
    values = brewer.pal(3, "Set2")[1:2] %>% 
      setNames(c("unk1", "unk2"))
  ) +
  labs(
    title = "CE chain length distribution",
    x = "Chain length",
    y = "Mole fraction phospholipids"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_acylcarbons.pdf"), width = 6, height = 6)
ggsave(here("05-png", "uw_acylcarbons.png"), width = 6, height = 6)


