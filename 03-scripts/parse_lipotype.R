library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))

# this directory should contain Lipid Maps spreadsheets as TSVs
dir_data   = here("01-rawdata", "lipotype")
# make sure the mole % file is picked up
pat_data   = "molp.csv" # filename pattern; UNIX glob
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T) %>% 
  # avoid MS autosave files
  .[which(!str_detect(., '~'))]

# how many datafiles did we get?
str_glue("Found {length(files_data)} datafiles:") %>% message()
files_data

# once parsed, data will be saved here for faster downstream loading
file_tidy  = here("02-tidydata", "lipotype_tidy.tsv")

# this block reads and binds multiple LipoType CSV files
# if you have just one file, that's fine, length(files_data) will =1.
# LipoType data files come with annotations baked in, so there is less to do here
# than for Lipid Maps and I use all standard functions, no custom helpers.
ltypdata_long = files_data %>% 
  lapply(., read_csv) %>% 
  bind_rows() %>%
  # there is a dud column and a couple dud rows that start with blank cells
  select(-fullname) %>% 
  filter(!is.na(feature)) %>% 
  # pivot to long format (eliminate sample-named columns)
  pivot_longer(cols = -names(cols_meta_ltyp), names_to = "eid", values_to = "frac_molar") %>%
  # make zeroes explicit
  mutate(frac_molar = as.numeric(frac_molar)) %>% 
  replace_na(list(frac_molar = 0)) %>% 
  # rename columns to match seelipids standard
  rename_at(names(cols_meta_ltyp), .funs = function(x){cols_meta_ltyp[x]}) %>% 
  # renormalize to mole fraction
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # remove zeroes to avoid plotting empty classes
  filter(frac_molar > 0) %>% 
  # rename some lipid classes to Lipid Maps standard
  mutate(
    class = ifelse(class == "DAG", "DG", class),
    class = ifelse(class == "TAG", "TG", class)
  ) %>% 
  # generate column for labeling plots
  separate(id, into = c(NA, "annot"), sep=' ', remove = FALSE)

# a little report
str_glue("Found {length(unique(ltypdata_long$id))} compounds in {length(unique(ltypdata_long$eid))} samples") %>% message()

# save parsed data to a new TSV
ltypdata_long %>% 
  write_tsv(file_tidy)

