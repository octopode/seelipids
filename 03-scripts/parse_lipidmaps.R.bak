library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))

# this directory should contain Lipid Maps spreadsheets as TSVs
dir_data   = here("01-rawdata", "lipidmaps")
pat_data   = "tsv" # filename pattern; UNIX glob
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T)

# how many datafiles did we get?
str_glue("Found {length(files_data)} datafiles:") %>% message()
files_data

file_tidy  = here("02-tidydata", "lipidmaps_tidy.tsv")

# Number of non-data lines to skip from the top of every raw TSV
lines_skip = 7

# collate all the Lipid Maps datafiles
lmapdata_wide = files_data %>%
  lapply(
    .,
    function(file_data){
      message(file_data)
      data_this_file <- file_data %>%
        # strip comment lines off the top
        read_tsv(skip = lines_skip) %>%
        # rename the extract ID column
        dplyr::rename(eid = `Sample ID`) %>%
        # replace "ND" with zero
        mutate_all(Vectorize(function(x){ifelse(x=="ND", 0, x)})) %>%
        # drop "invented" columns and rows
        select(-contains("...")) %>%
        # in these data files every column except the ID should be numeric
        mutate(across(!one_of("eid"), as.numeric)) %>% 
        drop_na()
      return(data_this_file)
    }
  ) %>%
  # more robust way to do the combined join/bind:
  bind_rows() %>% 
  group_by(eid) %>%
  summarise(across(everything(), function(x){c(na.omit(x), 0)[[1]]})) %>% 
  # strip the totals
  filter(!str_detect(eid, "Total"))

# pivot data to long format and parse the LM nomenclature
lmapdata_long = lmapdata_wide %>%
  # `rab` stands for "relative/raw abundance"
  # `id` always stands for *compound* id
  pivot_longer(cols = -one_of("eid"), names_to = "id", values_to = "rab") %>%
  # see *lipidomics_helpers.R to understand/check outputs
  parse_lipidmaps_id() %>%
  # normalize by total lipid abundance
  group_by(eid) %>%
  # na.rm is IMPORTANT because missing eid-id combinations are not explicit!
  mutate(frac_molar = rab/sum(rab, na.rm = TRUE)) %>%
  # dealing with compound non-overlap
  replace_na(list("frac_molar" = 0)) %>% 
  group_by(id)

# a little report
str_glue("Found {length(unique(lmapdata_long$id))} compounds in {length(unique(lmapdata_long$eid))} samples") %>% message()

# save parsed data to a new TSV
lmapdata_long %>% 
  write_tsv(file_tidy)
