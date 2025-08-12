
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

create_sourcemods = function(sample.sheet, id_col, group_col, gen.object) {
  
  # capture the column names as symbols
  id_col = rlang::ensym(id_col)
  group_col = rlang::ensym(group_col)
  
  ind.names = gen.object@ind.names
  
  sample.info.sub = sample.sheet %>%
    dplyr::filter(!!id_col %in% ind.names) %>%
    dplyr::slice(match(ind.names, .[[rlang::as_name(id_col)]])) %>%
    dplyr::select(!!id_col, !!group_col)
  
  return(sample.info.sub)
}

  
  