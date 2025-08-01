
library(magrittr)

#############################################################################################
# CREATE INDEX FILES
#############################################################################################

# (1) A data frame of the samples for each plate, where these two columns have to be present: 
# "sample_id" and "well" 

# Run this script separately for each plate

#############################################################################################

# Note that the i5 index has an extra AT
# and the i7 index has an extra T
# I manually added these, based on the sequencing results

#############################################################################################

create_internal_indexes = function(sample_ids){
  
  index_reference = data.frame(
    row = rep(LETTERS[1:8], times = 12),
    col = rep(1:12, each = 8),
    well = c(
      paste0(rep(LETTERS[1:8], times = 12), rep(1:12, each = 8))
    ),
    i5_index = rep(c(
      "CCGAATAT", "TTAGGCAAT", "AACTCGTCAT", "GGTCTACGTAT", "GATACCAT", "AGCGTTGAT", "CTGCAACTAT", "TCATGGTCAAT"
    ), times = 12),
    i7_index = rep(c(
      "CTAACGT", "TCGGTACT", "GATCGTTGT", "AGCTACACTT", "ACGCATT", "GTATGCAT", "CACATGTCT", "TGTGCACGAT", 
      "GCATCAT", "ATGCTGTT", "CATGACCTT", "TGCAGTGAGT"
    ), each = 8)
  )
  
  # generate a table that has allocated the correct internal index to each sample
  matched_samples = sample_ids %>%
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
