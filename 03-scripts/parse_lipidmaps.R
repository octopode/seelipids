library(tidyverse)
library(here)
library(readxl)

source(here("03-scripts", "seelipids_helpers.R"))

# this directory should contain Lipid Maps spreadsheets as TSVs
dir_data   = here("01-rawdata", "lipidmaps")
pat_data   = "xlsx" # filename pattern; UNIX glob
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T) %>% 
  # avoid MS autosave files
  .[which(!str_detect(., '~'))]

# how many datafiles did we get?
str_glue("Found {length(files_data)} datafiles:") %>% message()
files_data

# once parsed, data will be saved here for faster downstream loading
file_tidy  = here("02-tidydata", "lipidmaps_tidy.tsv")

# this block reads and binds multiple Lipid Maps Excel files
# if you have just one file, that's fine, length(files_data) will =1.
lmapdata_long = files_data %>% 
  # read_lmx() is the magic helper function here
  lapply(., read_lmx) %>% 
  bind_rows() %>% 
  select(-path, -sheet) %>% 
  # now parse and normalize
  # mark non-detection as zero
  replace_na(list(rab = 0)) %>% 
  # see seelipids_helpers.R to understand/check outputs
  parse_lipidmaps_id() %>% # needs to be sped up
  # normalize by total lipid abundance
  group_by(eid) %>%
  mutate(frac_molar = rab/sum(rab)) %>%
  # remove zeroes to avoid plotting empty classes
  filter(frac_molar > 0) %>% 
  group_by(id)

# a little report
str_glue("Found {length(unique(lmapdata_long$id))} compounds in {length(unique(lmapdata_long$eid))} samples") %>% message()

# save parsed data to a new TSV
lmapdata_long %>% 
  write_tsv(file_tidy)
