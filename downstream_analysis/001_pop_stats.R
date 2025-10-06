library(SNPfiltR)
library(vcfR)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(PopGenReport)
library(ggplot2)

popmap.broad = read.table("populations/broad_pops.txt", header = TRUE)
popmap.country = read.table("populations/country_pops.txt", header = TRUE)
popmap.site = read.table("populations/site_pops.txt", header = TRUE)

cryo.snps = vcfR::read.vcfR("populations/populations.snps.filtered_mac3.vcf")

# convert to genind
cryo.snps_genind = vcfR::vcfR2genind(cryo.snps) ;cryo.snps_genind
# genlight
cryo.snps_genlight = vcfR::vcfR2genlight(cryo.snps) ;cryo.snps_genlight@ind.names

dartR::gl2faststructure(x = cryo.snps_genlight, outpath=getwd(),
                        outfile = "populations/faststructure_input_mac3.str")


cryo.snps_genlight@pop<-as.factor(cryo.snps_genlight@ind.names)
sample.div = StAMPP::stamppNeisD(cryo.snps_genlight, pop = F)
StAMPP::stamppPhylip(distance.mat=sample.div, 
                     file="populations/splitstree_mac3.txt")

# Extract sample names from the genlight object
filtered_samples <- cryo.snps_genlight@ind.names
# Extract original sample names from popmap
original_samples <- popmap.broad$id  # assuming popmap.broad$id is a character vector
# Find missing samples
missing_samples <- setdiff(original_samples, filtered_samples)
# Print them
cat("Samples present in popmap but missing from genlight:\n")
print(missing_samples)

# Filter out those missing samples from each popmap
popmap.broad <- popmap.broad[!popmap.broad$id %in% missing_samples, ]
popmap.country <- popmap.country[!popmap.country$id %in% missing_samples, ]
popmap.site <- popmap.site[!popmap.site$id %in% missing_samples, ]

# create copies of the original data, but assign different population structures

# broad grouping
broad_grouping.genind = cryo.snps_genind
broad_grouping.genind@pop = as.factor(popmap.broad$pop)

broad_grouping.genlight = cryo.snps_genlight
broad_grouping.genlight@pop = as.factor(popmap.broad$pop)

# country grouping
country_grouping.genlight = cryo.snps_genlight
country_grouping.genlight@pop = as.factor(popmap.country$pop)

country_grouping.genind = cryo.snps_genind
country_grouping.genind@pop = as.factor(popmap.country$pop)

# site grouping
site_grouping.genind = cryo.snps_genind
site_grouping.genind@pop = as.factor(popmap.site$pop)

site_grouping.genlight = cryo.snps_genlight
site_grouping.genlight@pop = as.factor(popmap.site$pop)

###########################################################################
# FST VALUES
###########################################################################

inds_to_keep = which(!indNames(broad_grouping.genind) %in% c("HKSA1", "HKSA2", "HKSA3"))

# remove those three HKSA samples 
broad_grouping.genind = broad_grouping.genind[inds_to_keep, ]
country_grouping.genind = country_grouping.genind[inds_to_keep, ]

# BROAD

broad.FST = hierfstat::genet.dist(broad_grouping.genind, method = "WC84")
broad.FST = as.data.frame(as.matrix(broad.FST)) %>%
  round(., digits = 2)

##############################################################################
# POP STATS

se <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}
##############################################################################

# BROAD
##############################################################################

broad_poppr_result = poppr::poppr(broad_grouping.genind)
broad_poppr_result <- broad_poppr_result[-nrow(broad_poppr_result), ]

broad.stats = hierfstat::basic.stats(data = broad_grouping.genind)
broad.AR = hierfstat::allelic.richness(data = broad_grouping.genind)

broad_poppr_result$Fis = colMeans(broad.stats$Fis, na.rm = TRUE)
broad_poppr_result$Fis_se   <- apply(broad.stats$Fis, 2, se)

broad_poppr_result$Ho = colMeans(broad.stats$Ho, na.rm = TRUE)
broad_poppr_result$Ho_se   <- apply(broad.stats$Ho, 2, se) 

broad_poppr_result$Hs = colMeans(broad.stats$Hs, na.rm = TRUE)
broad_poppr_result$Hs_se   <- apply(broad.stats$Hs, 2, se) 

broad_poppr_result$allelerichness = colMeans(broad.AR$Ar, na.rm = TRUE)
broad_poppr_result$allelerichness_se = apply(broad.AR$Ar, 2, se)

broad_private_alleles = poppr::private_alleles(broad_grouping.genind, report = "data.frame",
                                               level = "population", count.alleles = FALSE) 

broad_total_private_perpop = aggregate(count ~ population, 
                                       data = broad_private_alleles, sum) 

broad_poppr_result$private_alleles = broad_total_private_perpop$count

##############################################################################
# COUNTRY
##############################################################################
country_poppr_result = poppr::poppr(country_grouping.genind)
country_poppr_result <- country_poppr_result[-nrow(country_poppr_result), ]

country.stats = hierfstat::basic.stats(data = country_grouping.genind)
country.AR = hierfstat::allelic.richness(data = country_grouping.genind)

country_poppr_result$Fis = colMeans(country.stats$Fis, na.rm = TRUE)
country_poppr_result$Fis_se   <- apply(country.stats$Fis, 2, se)

country_poppr_result$Ho = colMeans(country.stats$Ho, na.rm = TRUE)
country_poppr_result$Ho_se   <- apply(country.stats$Ho, 2, se) 

country_poppr_result$Hs = colMeans(country.stats$Hs, na.rm = TRUE)
country_poppr_result$Hs_se   <- apply(country.stats$Hs, 2, se) 

country_poppr_result$allelerichness = colMeans(country.AR$Ar, na.rm = TRUE)
country_poppr_result$allelerichness_se = apply(country.AR$Ar, 2, se)

country_private_alleles = poppr::private_alleles(country_grouping.genind, report = "data.frame",
                                                 level = "population", count.alleles = FALSE) 

country_total_private_perpop = aggregate(count ~ population, 
                                         data = country_private_alleles, sum)

country_poppr_result$private_alleles = country_total_private_perpop$count

##############################################################################
# SITE
##############################################################################

site_poppr_result = poppr::poppr(site_grouping.genind)
site_poppr_result <- site_poppr_result[-nrow(site_poppr_result), ]

site.stats = hierfstat::basic.stats(data = site_grouping.genind)
site.AR = hierfstat::allelic.richness(data = site_grouping.genind)

site_poppr_result$Fis = colMeans(site.stats$Fis, na.rm = TRUE)
site_poppr_result$Fis_se   <- apply(site.stats$Fis, 2, se)

site_poppr_result$Ho = colMeans(site.stats$Ho, na.rm = TRUE)
site_poppr_result$Ho_se   <- apply(site.stats$Ho, 2, se) 

site_poppr_result$Hs = colMeans(site.stats$Hs, na.rm = TRUE)
site_poppr_result$Hs_se   <- apply(site.stats$Hs, 2, se) 

site_poppr_result$allelerichness = colMeans(site.AR$Ar, na.rm = TRUE)
site_poppr_result$allelerichness_se = apply(site.AR$Ar, 2, se)

site_private_alleles = poppr::private_alleles(site_grouping.genind, report = "data.frame",
                                              level = "population", count.alleles = FALSE) 

site_total_private_perpop = aggregate(count ~ population, 
                                      data = site_private_alleles, sum)

site_poppr_result$private_alleles = site_total_private_perpop$count

##############################################################################
# WRITE OUT FILES
##############################################################################

write.csv(broad_poppr_result, 
          "populations/broad_poppr_results.csv", row.names = FALSE)
write.csv(country_poppr_result, 
          "populations/country_poppr_results.csv", row.names = FALSE)
write.csv(site_poppr_result, 
          "populations/site_poppr_results.csv", row.names = FALSE)

##############################################################################
# DARTR experimenting 
##############################################################################

cryo.snps = vcfR::read.vcfR("populations/populations.snps.filtered_mac3.vcf")
cryo.snps_genlight = vcfR::vcfR2genlight(cryo.snps)

pcoa = dartR::gl.pcoa(broad_grouping.genlight, plot.out = F, verbose = 0)
pcoaplot = dartR::gl.pcoa.plot(pcoa, 
                         broad_grouping.genlight, 
                         verbose = 0)

# run Structure

cryst.struc = dartR::gl.run.structure(broad_grouping.genlight,
                                      k.range = 2:6,
                                      num.k.rep = 10, 
    exec = "D:/structure_windows_console/console/structure.exe",
                                      plot.out = FALSE,
                                      burnin = 50000,
                                      numreps = 100000,
                                      noadmix = FALSE )

saveRDS(cryst.struc, "path.rds")
