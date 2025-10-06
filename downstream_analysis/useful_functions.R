
library(magrittr)

#############################################################################################
# CREATE INDEX FILES
#############################################################################################

# (1) A data frame of the samples for each plate, where these two columns have to be present: 
# "sample_id" and "well" 

# Run this script separately for each plate

#############################################################################################
# Internal indexes are hardcoded into this script, using:
# ClaI and EcoRI iTru i5 and iTru i7 from Bayona-Vasquez et al., 2019
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6791345/
# Adapterama III: Quadruple-indexed, double/triple-enzyme RADseq libraries (2RAD/3RAD)

# Note that the i5 index has an extra AT
# and the i7 index has an extra T

# I manually added these, based on the sequencing results

#############################################################################################

create_internal_indexes = function(plate_samples){
  
  index_reference = data.frame(
    row = rep(LETTERS[1:8], times = 12),
    col = rep(1:12, each = 8),
    well = c(
      paste0(rep(LETTERS[1:8], times = 12), rep(1:12, each = 8))
    ),
	# Rows A -> H on a 96-well plate
    i5_index = rep(c(
      "CCGAATAT", "TTAGGCAAT", "AACTCGTCAT", "GGTCTACGTAT", "GATACCAT", "AGCGTTGAT", "CTGCAACTAT", "TCATGGTCAAT"
    ), times = 12),
	# Columns 1 -> 12 on a 96-well plate
    i7_index = rep(c(
      "CTAACGT", "TCGGTACT", "GATCGTTGT", "AGCTACACTT", "ACGCATT", "GTATGCAT", "CACATGTCT", "TGTGCACGAT", 
      "GCATCAT", "ATGCTGTT", "CATGACCTT", "TGCAGTGAGT"
    ), each = 8)
  )
  
  # generate a table that has allocated the correct internal index to each sample
  matched_samples = plate_samples %>%
    dplyr::left_join(index_reference, by = "well") %>%
    dplyr::select(i5_index, i7_index, sample_id)
  
  return(matched_samples)
  
}

# write out result as txt

# write.table(matched_samples, 
#            file = "sample_info/dean_guerich_included/internal_indexes_plate_2.txt", 
#            sep = "\t",             # Use a tab space as the separator
#            row.names = FALSE,     # Do not write row names
#            col.names = FALSE,      # Write column names
#            quote = FALSE,         # Avoid quoting character data
#            eol = "\n")            # Ensure Unix line endings

#############################################################################################
# CREATE POP FILE
#############################################################################################

# INPUTS

# (1) A data frame of all your samples, where these two columns have to be present: 
# "sample_id" and "origin" 

create_pop_file = function(sample_sheet){
  
  pop_info = sample_sheet %>%
    janitor::clean_names()  %>%
    dplyr::select(sample_id, origin)
  
  return(pop_info)
}

# write.table(pop_info, 
#             file = "pops_all.txt", 
#             sep = "\t",             # Use a tab space as the separator
#             row.names = FALSE,     # Do not write row names
#             col.names = FALSE,      # Write column names
#             quote = FALSE,         # Avoid quoting character data
#             eol = "\n")            # Ensure Unix line endings

#############################################################################################
# CREATE SOURCE MODIFIERS FILE - to accompany a genlight or genind object
# this makes sure that all the samples in the genlight/genind have an associated file with
# population info, as some samples are dropped during filtering etc., and the tallies won't
# always match up

# useful for downstream processing, and uploading to the SRA database later
#############################################################################################

# INPUTS

# (1) From your sample sheet: the column of sample IDs and the column for the grouping variable of interest
# (2) The genlight or genind object containing your SNP data

create_sourcemods = function(sample.sheet, id_col, gen.object, group_col = NULL, 
                              keep_all = FALSE, lat_col = NULL, lon_col = NULL) {
  
  # Convert strings to symbols
  id_col = rlang::sym(id_col)
  if (!is.null(group_col)) group_col = rlang::sym(group_col)
  if (!is.null(lat_col)) lat_col = rlang::sym(lat_col)
  if (!is.null(lon_col)) lon_col = rlang::sym(lon_col)
  
  ind.names = gen.object@ind.names
  
  sample.info.sub = sample.sheet %>%
    dplyr::filter(!!id_col %in% ind.names) %>%
    dplyr::slice(match(ind.names, .[[rlang::as_string(id_col)]]))
  
  if (!keep_all) {
    if (is.null(group_col)) stop("If keep_all = FALSE, you must provide a group_col.")
    sample.info.sub = sample.info.sub %>%
      dplyr::select(!!id_col, !!group_col)
  } else {
    if (!is.null(lat_col) && !is.null(lon_col)) {
      lat_name = rlang::as_string(lat_col)
      lon_name = rlang::as_string(lon_col)
      
      if (all(c(lat_name, lon_name) %in% names(sample.info.sub))) {
        sample.info.sub = sample.info.sub %>%
          dplyr::mutate(
            # Convert to numeric safely here
            lat_num = as.numeric(.data[[lat_name]]),
            lon_num = as.numeric(.data[[lon_name]])
          ) %>%
          dplyr::mutate(
            coordinates = paste0(
              formatC(abs(lat_num), format = "f", digits = 2), " ",
              ifelse(lat_num >= 0, "N", "S"), " ",
              formatC(abs(lon_num), format = "f", digits = 2), " ",
              ifelse(lon_num >= 0, "E", "W")
            )
          ) %>%
          dplyr::select(-lat_num, -lon_num)  # remove helper cols
      } else {
        warning("Latitude or longitude columns not found; 'coordinates' column not created.")
      }
    }
  }
  
  return(sample.info.sub)
}



  
################################################################
# PLOT STRUCTURE/fastSTRUCTURE output
################################################################

# inputs:
# output file from fastSTRUCTURE (.meanQ)
# sample info data frame, with two columns: id and pop

#################################################################

structure.plot = function(strucfile, sampleinfo, kval="N"){
  
  f = readLines(strucfile, warn = FALSE)
  
  qmat = read.delim(
    strucfile, 
    header = FALSE,
    sep = "",
    strip.white = TRUE
  )
  
  names(qmat) = c(paste0("pop_",1:(ncol(qmat))))
  
  qmat$name = sampleinfo$id
  qmat$population = sampleinfo$pop
  
  qmat$population = as.factor(qmat$population)
  qmat = qmat %>% 
    pivot_longer(c(-name, -population), names_to = "group", values_to = "probability")
  
  #my_pal = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  my_pal = c("black", "#D95F02", "#1B9E77", "yellow","#666666", "#7570B3", "#66A61E", "#A6761D",
             "#E7298A", "red", "lightblue", "royalblue")
  
  struc.plot = qmat %>% 
    ggplot(aes(x = name, y = probability, fill = group)) +
    geom_bar(stat = "identity", width = 1.0) +
    theme_bw() +
    labs(title = paste0("K = ", kval)) +
    ylab("Probability") +
    xlab("") +
    #guides(fill=guide_legend(title="Membership")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off") +
    scale_fill_manual(values=c(my_pal)) +
    facet_grid(~population, scales = "free_x", space = "free" )
  
  return(struc.plot)
  
}#structure.plot



hapmap.global = function(struc.file, sample.info, gps.data, piesize = 0.85){
  
  f = readLines(struc.file, warn = FALSE)
  
  qmat = read.delim(
    struc.file, 
    header = FALSE,
    sep = "",
    strip.white = TRUE
  )
  
  names(qmat) = c(paste0("pop_",1:(ncol(qmat))))
  
  qmat$name = sample.info$id
  qmat$population = sample.info$pop
  
  qmat$population = as.factor(qmat$population)
  qmat = qmat %>% 
    pivot_longer(c(-name, -population), names_to = "group", values_to = "probability")
  
  #my_pal = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  my_pal = c("black", "#D95F02", "#1B9E77", "yellow","#666666", "#7570B3", "#66A61E", "#A6761D",
             "#E7298A", "red")
  
  qmat_merged = qmat %>%
    group_by(name) %>%
    mutate(lat = gps.data$lat[match(name, gps.data$name)],
           lon = gps.data$lon[match(name, gps.data$name)]) %>%
    ungroup()
  
  qmat_merged$lat = as.numeric(qmat_merged$lat)
  qmat_merged$lon = as.numeric(qmat_merged$lon)
  
  # this is a summary for each individual sample
  qmat_wide = qmat_merged %>%
    tidyr::pivot_wider(names_from = group, values_from = probability)
  
  # this averages all the data per unique group, across all samples for that group
  qmat_averaged = qmat_wide %>%
    group_by(population) %>%
    summarize(
      lat = first(lat),      # Keep the first value of lat
      lon = first(lon),      # Keep the first value of lon
      across(starts_with("pop"), mean, na.rm = TRUE)  # Average pop_* columns
    ) %>%
    dplyr::rename(name = population) 
  
  qmat_averaged = qmat_averaged %>%
    mutate(
      lat = ifelse(is.na(lat), 35.0, lat),  # Replace NA lat with Cyprus lat
      lon = ifelse(is.na(lon), 33.0, lon)   # Replace NA lon with Cyprus lon
    )
  
  world_map = rnaturalearth::ne_countries(
    scale = "medium", 
    returnclass = "sf"
  ) 
  
  pies = ggplot() +
    
    geom_sf(data = world_map, alpha = 0.2, fill = "grey80") +
    
    scatterpie::geom_scatterpie(
      #data = na.omit(qmat_wide), # plot pies for each individual sample
      data = na.omit(qmat_averaged), # plot pies representing the average per group
      aes(x = lon, y = lat, group = name),
      #aes(x = lon, y = lat, group = name),
      cols = pop_cols = grep("^pop_", colnames(qmat_wide), value = TRUE),
      #alpha = 0.3,
      #color = NA, # Remove borders for pie charts
      color = "black", # Remove borders for pie charts
      linewidth = 0.2,
      pie_scale = piesize # Adjust pie size if necessary
    ) +
    
    # uncomment if you want labels on the pies
    # geom_text(
    #   data = na.omit(qmat_averaged),
    #   aes(x = lon, y = lat, label = name), # Assuming 'name' is the group label
    #   color = "black",
    #   fontface = "bold",
    #   nudge_y = 0.1 # Adjust this value if needed
    # ) +
    
    scale_fill_manual(values = my_pal) +
    
    coord_sf(
      ylim = c(95,-60),
      crs = 4326,
      expand = FALSE
    ) + 
    
    xlab("Longitude") + 
    ylab("Latitude") +
    theme_classic() +
    theme(legend.position = "none")
  
  return(pies)
} #hapmap.global
