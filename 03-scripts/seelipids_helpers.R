# HELPER FUNCTIONS for lipidomics analyses and plotting

library(here)
library(tidyverse)
library(readxl)
library(ggpubr)
library(RColorBrewer)

## IMPORTANT MAPPINGS

# All the lipid classes currently recognized by seelipids.
# These are all acceptable in input, but some are synonymous and will be
# mapped to a more standardized nomenclature.
classes = c(
  "P-PC", "O-PC", "PC", "LPC",
  "P-PE", "O-PE", "PE", "LPE",
  "PS", "LPS",
  "PI", "LPI",
  "CL", "PG", "LPG",
  "PA",
  "Cer", "HexCer", "HexCER", "LacCer", "LacCER",  "Cer1P", "SM",
  "FA", "DG", "TG",
  "AC", "CE"
)

## LIPID MAPS PARSING HELPERS

# Read multi-sheet Lipid Maps Excel file,
# return raw abundances in long format.
# Compatible with LM format as of 20220409.
read_lmx = function(files, skip = 9, na = "ND"){
  # get the data sheets in the files
  tibble(path = files) %>% 
    rowwise() %>% 
    mutate(sheet = path %>% excel_sheets() %>% list()) %>% 
    unnest(sheet) %>% 
    filter(str_detect(sheet, "(detected|data)")) %>% # backward-compatible w/2020 format
    rowwise() %>% 
    mutate(
      # read all the data out
      data = read_excel(
        path, 
        sheet,
        na = na,
        skip = skip
      ) %>% 
        select(!contains("...")) %>% 
        rename(eid = `Sample ID`) %>%
        list()
    ) %>% 
    unnest(data) %>% 
    # get rid of "Total" rows
    filter(
      !is.na(eid) &
        !str_detect(eid, "Total")
    ) %>% 
    # coerce to ensure meltability
    # keeps path and sheet metadata, which can come in handy
    mutate(across(-c(eid, path, sheet), as.numeric)) %>% 
    pivot_longer(
      cols = -c(eid, path, sheet),
      names_to  = "id",
      values_to = "rab"
    )
}

# parse Lipid Species CSV file of the format sent by ETS @ UW
# Can read multiple files; each sample is expected to occur _once_
# in the entire file set
read_uw = function(files){
  # load all files and add a column with source filename
  files %>% 
    read_csv(id = "fname") %>% 
    rename(eid = `Sample ID`) %>% 
    # pivot to long format
    # if one EID occurs in multiple files, the entries can be
    # distinguished using the fname column
    group_by(fname, eid) %>%
    pivot_longer(
      cols = -c(fname, eid), 
      names_to = "id", 
      # rab stands for "relative abundance"
      values_to = "rab"
    )
}

# parse Lipid Species Excel spreadsheet of the format sent by ETS @ UW
# Can read multiple files; each sample is expected to occur _once_
# in the entire file set
read_uwx = function(files){
  # load all files and add a column with source filename
  lapply(
    files,
    function(x){read_excel(x) %>% mutate(fname = basename(x))}
  ) %>% 
    # bind the files together
    as.list() %>% 
    bind_rows() %>% 
    rename(eid = `Sample ID`) %>% 
    # pivot to long format
    # if one EID occurs in multiple files, the entries can be
    # distinguished using the fname column
    group_by(fname, eid) %>%
    pivot_longer(
      cols = -c(fname, eid), 
      names_to = "id", 
      # rab stands for "relative abundance"
      values_to = "rab"
    )
}

# just the standard error of the mean
serr = function(x){sd(x)/sqrt(length(x))}

# How many acyl chains does passed class of lipids have?
class2tails = Vectorize(function(class){
  if(str_detect(class, 'L|AC|CE|FA')){
    return(1L)
  }
  if(str_detect(class, 'T')){
    return(3L)
  }else{
    return(2L)
  }
})

# Take a long-format tibble with an `id` column containing lipid maps compound names,
# break that column out into several informative columns
#NTS 20220502: May need to be updated for changes in nomenclature between 2020, 2022.
#NTS 20220909: This is kludgy and could stand to be refactored with stringr...but it works!
parse_lipidmaps_id = function(longdf){
  longdf %>% 
    group_by(id) %>%
    mutate(
      # parse headgroup-level structural data
      struc  = id %>%
        str_split('[ _:/]'),
      # headgroup class
      class =  unlist(struc)[[1]] %>%
        # this refers to the shared list of recognized headgroups
        factor(levels = classes),
      # tails
      tails = class2tails(class),
      # total # hydroxy groups
      hydrox = unlist(struc)[[2]] %>%
        str_extract("[nmdt]") %>%
        c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
      # total carbons
      carbon = unlist(struc)[[2]] %>%
        str_remove("[nmdt]") %>%
        as.integer(),
      # total double bonds
      dbonds  = unlist(struc)[[3]] %>%
        as.integer(),
      rt      = unlist(struc)[[4]] %>%
        as.numeric(),
      # fatty acid-level structural data
      hydrsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[6]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[6]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[7]] %>%
          as.integer(),
        NA
      ),
      hydrsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[8]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[8]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[9]] %>%
          as.integer(),
        NA
      ),
      hydrsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[10]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[10]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[11]] %>%
          as.integer(),
        NA
      ),
      # sn1/2 positions resolved?
      stereo = str_detect(id, '/') %>% unlist(),
      # make an annotation string that is as specific as possible
      annot = ifelse(
        '|' %in% unlist(struc),
        # if there is acyl chain info
        id %>%
          str_split(' \\| ') %>%
          unlist() %>%
          .[[2]],
        # if there's no acyl chain info
        paste(unlist(struc)[[2]], unlist(struc)[[3]], sep=':')
      )
    )
}

## LIPOTYPE PARSING HELPERS

# non-sample columns in a lipotype datafile
# also map to their names in my own convention
cols_meta_ltyp = c(
  "feature" = "id", 
  "class" = "class", 
  "structuralcategory" = "structuralcategory", 
  "functionalcategory" = "functionalcategory", 
  "totallength" = "carbon", 
  "totaldb" = "dbonds", 
  "totaloh" = "hydrox"
)

# UW format parsing helpers
parse_uw_id = function(longdf){
  longdf %>% 
    group_by(id) %>%
    mutate(
      # use separate()? 20250303
      # parse headgroup-level structural data
      struc  = id %>%
        str_split('[ _:/]'),
      # headgroup class
      class =  unlist(struc)[[1]] %>%
        # this refers to the shared list of recognized headgroups
        factor(levels = classes),
      # tails
      tails = class2tails(class),
      # total # hydroxy groups
      hydrox = unlist(struc)[[2]] %>%
        str_extract("[nmdt]") %>%
        c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
      # total carbons
      carbon = unlist(struc)[[2]] %>%
        str_remove("[nmdt]") %>%
        as.integer(),
      # total double bonds
      dbonds  = unlist(struc)[[3]] %>%
        as.integer(),
      rt      = unlist(struc)[[4]] %>%
        as.numeric(),
      # fatty acid-level structural data
      hydrsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[6]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[6]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn1 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[7]] %>%
          as.integer(),
        NA
      ),
      hydrsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[8]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[8]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn2 = ifelse(
        "|" %in% unlist(struc),
        unlist(struc)[[9]] %>%
          as.integer(),
        NA
      ),
      hydrsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[10]] %>%
          str_extract("[nmdt]") %>%
          c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
        NA
      ),
      carbsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[10]] %>%
          str_remove("[nmdt]") %>%
          as.numeric(),
        NA
      ),
      dbonsn3 = ifelse(
        ("|" %in% unlist(struc)) & (tails == 3),
        unlist(struc)[[11]] %>%
          as.integer(),
        NA
      ),
      # sn1/2 positions resolved?
      stereo = str_detect(id, '/') %>% unlist(),
      # make an annotation string that is as specific as possible
      annot = ifelse(
        '|' %in% unlist(struc),
        # if there is acyl chain info
        id %>%
          str_split(' \\| ') %>%
          unlist() %>%
          .[[2]],
        # if there's no acyl chain info
        paste(unlist(struc)[[2]], unlist(struc)[[3]], sep=':')
      )
    )
}

## PLOTTING HELPERS

# full 14-headgroup color map
# 20220919 REVISED for Lipid Maps/LipoType cross-compatibility
chroma_cl = c(
  brewer.pal(4, "Blues")[2:4], #PCs
  brewer.pal(4, "Greens")[2:4], #PEs
  brewer.pal(6, "YlOrRd")[3:2], #PAs; beware that these are very similar to DG color...
  brewer.pal(12, "Paired")[5:12] #PSs, DG, PI...
) %>%
  setNames(
    c(
      "P-PC", "PC", "LPC",
      "P-PE", "PE", "LPE",
      "PA", "LPA",
      "PS", "LPS",
      "DG",
      "PI",
      "SM",
      "Cer",
      "PG", "TG"
    )
  ) %>% 
  # Some additional classes with manual color spec
  c(., c(
    LPI   = brewer.pal(9, "Oranges")[[7]],
    CL    = "darkred",
    Cer1P = "mediumorchid",
    AC    = "grey50",
    CE    = "grey25",
    EE    = "grey25"
  )) %>% 
  # The order of the named list matters for legend construction.
  # CLASSES NOT INCLUDED HERE WILL NOT BE PLOTTED!
  # This is the order in which bars will be stacked;
  # it is loosely based on intrinsic curvature:
  .[c("P-PC", "PC", "LPC", "P-PE", "PE", "LPE", "PS", "LPS", "CL", "PI", "LPI",
      "PA", "LPA", "PG", "Cer", "Cer1P", "SM", "DG", "TG", "AC", "CE", "EE")]

# colors for label text on a bar plot
chroma_tx = c(
  "PC"="black",  "P-PC"="black", "LPC"="white",
  "PE"="black",  "P-PE"="black", "LPE"="white",
  "PS"="black",  "LPS"="white",
  "DG"="black",
  "PI"="black",  "LPI"="black",
  "PA"="black",  "LPA"="black",
  "CL"="white",
  "SM"="black",
  "Cer"="white", "Cer1P"="white",
  "PG"="black",  "TG"="black",
  "AC"="white",  "CE"="white", "EE"="white"
)

# break generator for acyl chain lengths
even_breaks = function(x, n=5){
  x = sort(round(x))
  if((x[[1]] %% 2) > 0){x[[1]] <- x[[1]]+1}
  spacing <- 2
  breaks <- seq(min(x), max(x), spacing)
  while(length(breaks) > n){
    spacing <- spacing + 2
    breaks <- seq(min(x), max(x), spacing)
  }
  return(breaks)
}

# WonB theme for slide figures
theme_pubk = function(...){
  theme_pubr(...) +
    theme(
      # axis options
      axis.line  = element_line(color = "white"),
      axis.ticks = element_line(color = "white"),
      axis.text  = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      # legend
      legend.background = element_blank(),
      legend.key   = element_rect(color = NA,  fill = NA),
      legend.text  = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      # panel
      panel.background = element_blank(),
      #panel.border = element_blank()
      # facetting
      strip.background = element_rect(fill = NA, color = "white"),
      strip.text       = element_text(color = "white"),
      # plot options
      plot.background = element_blank(),
      plot.title = element_text(color = "white")
    )
}

# custom stacked-barplot wrapper for headgroups
gg_headgp = function(
  data,
  darkmode = FALSE,
  thres_draw = 0, # can pass a mole fraction threshold below which color blocks are removed.
  label_frac = 0.015, # min fraction at which bar get labeled; set to 1 to disable labels.
  label_size = 1.5,
  # aesthetics get passed in as naked args
  ...
){
  # unpack the ellipsis args as strings
  mapstrs = lapply(rlang::enexprs(...), as.character)
  # ensure ordering
  data = data %>% arrange(class)
  this_gg = data %>%
    ## apply threshold
    filter(eval(sym(mapstrs$y)) >= thres_draw) %>% 
    ## EXAMPLE softcoded to use whatever is passed as x
    #summarize(wavg = sum(eval(sym(mapstrs$x))*frac_molar)/sum(frac_molar)),
    ## and renorm
    group_by(eval(sym(mapstrs$x))) %>% 
    mutate(across(mapstrs$y, function(x){x/sum(x)})) %>% 
    # cut the factor down to include only the classes in the data
    mutate(class = class %>% factor(., levels = unique(as.character(.)))) %>% 
    ggplot(mapping = aes(...)) +
    # draw little hairlines in between compounds
    geom_col(
      width = 0.95, 
      size = 0.05, 
      color = ifelse(darkmode, "black", "white")
    )
  #print(this_gg)
  # label individual compounds if they are resolved in data
  if("annot" %in% colnames(data)){
    this_gg = this_gg +
      geom_text(
        aes(
          label = ifelse(frac_molar >= label_frac, annot, ''),
          color = class
        ),
        size = label_size,
        position = position_stack(vjust=0.5)
      )
  }
  if(darkmode){
    this_gg = this_gg + theme_pubk()
  }else{
    this_gg = this_gg + theme_pubr()
  }
  # finishing touches
  this_gg + 
    # assign scales
    # strip out unused classes from the scales
    scale_fill_manual( values = chroma_cl[data$class %>% as.character() %>% unique()]) +
    scale_color_manual(values = chroma_tx[data$class %>% as.character() %>% unique()]) +
    # position the legend and text
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle=90, vjust=0.5)
    ) +
    guides(color = "none")
}

# custom barplot array wrapper for acyl chains
gg_acylch = function(
  data,
  darkmode = FALSE,
  meanline = TRUE,
  alpha_odd = 1, # change to give odd x-coords a lighter alpha
  # can override default color mapping
  chroma = chroma_cl,
  # facet_grid() specification
  rows = vars(sp),
  cols = vars(class),
  # aesthetics get passed in as naked args
  ...
){
  # unpack the ellipsis args as strings
  mapstrs = lapply(rlang::enexprs(...), as.character)
  # start building plot
  this_gg = data %>% 
    ggplot(aes(...)) +
    facet_grid(
      rows = rows, 
      cols = cols
    ) +
    geom_col(
      # handles shading of odd chains if desired
      aes(alpha = ifelse(stat(x) %% 2, alpha_odd, 1)),
      size = 0.05, 
      color = ifelse(darkmode, "black", "white")
    )
  if(meanline){
    this_gg = this_gg +
      geom_vline(
        data = data %>% 
          # calc weighted average
          # workaround to play nice with vars() arg
          group_by_(
            str_sub(as.character(rows), 2L), 
            str_sub(as.character(cols), 2L)
          ) %>% 
          # softcoded to use whatever is passed as x
          summarize(wavg = sum(eval(sym(mapstrs$x))*frac_molar)/sum(frac_molar)),
        aes(xintercept = wavg),
        color = ifelse(darkmode, "white", "black"),
        size = 1,
        linetype="dashed"
      )
  }
  if(darkmode){
    this_gg = this_gg + theme_pubk()
  }else{
    this_gg = this_gg + theme_pubr()
  }
  # finishing touches
  this_gg +
    scale_x_continuous(breaks = even_breaks) +
    scale_y_continuous() +
    scale_fill_manual(values = chroma) +
    #scale_alpha(range = c(0.4, 2)) +
    scale_alpha_identity() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    guides(alpha="none", fill="none")
}
