library(SNPfiltR)
library(vcfR)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(PopGenReport)
library(ggplot2)

popmap.broad = read.table("populations/broad_pops.txt", header = TRUE)
popmap.broad

popmap.country = read.table("populations/country_pops.txt", header = TRUE)
popmap.country

popmap.site = read.table("populations/site_pops.txt", header = TRUE)
popmap.site

cryo.snps = vcfR::read.vcfR("populations/populations.snps.filtered.vcf")

# convert to genind
cryo.snps_genind = vcfR::vcfR2genind(cryo.snps)
cryo.snps_genind

# convert to genlight
cryo.snps_genlight = vcfR::vcfR2genlight(cryo.snps)
cryo.snps_genlight@ind.names



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
# these were the low quality samples filtered out earlier
popmap.broad <- popmap.broad[!popmap.broad$id %in% missing_samples, ]
popmap.country <- popmap.country[!popmap.country$id %in% missing_samples, ]
popmap.site <- popmap.site[!popmap.site$id %in% missing_samples, ]

write.table(popmap.broad, file = "populations/broad_pops_filtered.txt", quote= F, row.names = F)
write.table(popmap.country, file = "populations/country_pops_filtered.txt", quote= F, row.names = F)
write.table(popmap.site, file = "populations/site_pops_filtered.txt", quote= F, row.names = F)

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

##################################
# population statistics
##################################

##########################
# get broad pop stats
##########################

broad_poppr_result = poppr::poppr(broad_grouping.genind)
broad_poppr_result <- broad_poppr_result[-nrow(broad_poppr_result), ] # remove last "total" row

country_poppr_result = poppr::poppr(country_grouping.genind)
country_poppr_result <- country_poppr_result[-nrow(country_poppr_result), ]

site_poppr_result = poppr::poppr(site_grouping.genind)
site_poppr_result <- site_poppr_result[-nrow(site_poppr_result), ]

#############
# add more stats to the above DFs

broad.stats = hierfstat::basic.stats(data = broad_grouping.genind)
country.stats = hierfstat::basic.stats(data = country_grouping.genind)
site.stats = hierfstat::basic.stats(data = site_grouping.genind)

broad.AR = hierfstat::allelic.richness(data = broad_grouping.genind)
country.AR = hierfstat::allelic.richness(data = country_grouping.genind)
site.AR = hierfstat::allelic.richness(data = site_grouping.genind)

# add to existing DFs
broad_poppr_result$Fis = colMeans(broad.stats$Fis, na.rm = TRUE)
broad_poppr_result$Ho = colMeans(broad.stats$Ho, na.rm = TRUE)
broad_poppr_result$Hs =colMeans(broad.stats$Hs, na.rm = TRUE) 

broad_poppr_result$allelerichness = colMeans(broad.AR$Ar, na.rm = TRUE)

country_poppr_result$Fis = colMeans(country.stats$Fis, na.rm = TRUE)
country_poppr_result$Ho = colMeans(country.stats$Ho, na.rm = TRUE)
country_poppr_result$Hs = colMeans(country.stats$Hs, na.rm = TRUE)

country_poppr_result$allelerichness = colMeans(country.AR$Ar, na.rm = TRUE)

site_poppr_result$Fis = colMeans(site.stats$Fis, na.rm = TRUE)
site_poppr_result$Ho = colMeans(site.stats$Ho, na.rm = TRUE)
site_poppr_result$Hs = colMeans(site.stats$Hs, na.rm = TRUE)

site_poppr_result$allelerichness = colMeans(site.AR$Ar, na.rm = TRUE)

# find private alleles and add to existing DFs

broad_private_alleles = poppr::private_alleles(broad_grouping.genind, report = "data.frame",
                                               level = "population", count.alleles = FALSE) 

broad_total_private_perpop = aggregate(count ~ population, 
                                       data = broad_private_alleles, sum) 

country_private_alleles = poppr::private_alleles(country_grouping.genind, report = "data.frame",
                                                 level = "population", count.alleles = FALSE) 

country_total_private_perpop = aggregate(count ~ population, 
                                         data = country_private_alleles, sum)

site_private_alleles = poppr::private_alleles(site_grouping.genind, report = "data.frame",
                                              level = "population", count.alleles = FALSE) 

site_total_private_perpop = aggregate(count ~ population, 
                                      data = site_private_alleles, sum)

broad_poppr_result$private_alleles = broad_total_private_perpop$count
country_poppr_result$private_alleles = country_total_private_perpop$count
site_poppr_result$private_alleles = site_total_private_perpop$count

write.csv(broad_poppr_result, 
          "populations/broad_poppr_results.csv", row.names = FALSE)
write.csv(country_poppr_result, 
          "populations/country_poppr_results.csv", row.names = FALSE)
write.csv(site_poppr_result, 
          "populations/site_poppr_results.csv", row.names = FALSE)


# FST values
broad.FST = hierfstat::genet.dist(broad_grouping.genind, method = "WC84")
broad.FST = as.data.frame(as.matrix(broad.FST))
write.csv(broad.FST, "populations/broadFST.csv", row.names = T)

country.FST = hierfstat::genet.dist(country_grouping.genind, method = "WC84")
country.FST = as.data.frame(as.matrix(country.FST))
write.csv(country.FST, "populations/countryFST.csv", row.names = T)

site.FST = hierfstat::genet.dist(site_grouping.genind, method = "WC84")
site.FST = as.data.frame(as.matrix(site.FST))
write.csv(site.FST, "populations/siteFST.csv", row.names = T)

#################################
# plotting
#################################

broad_poppr_result = read.csv("populations/broad_poppr_results.csv") %>%
  dplyr::select(-c(File, Hexp, SE, E.5, eMLG)) %>%
  dplyr::rename("PA" = "private_alleles",
                "AR" = "allelerichness") %>%
  tidyr::pivot_longer(cols = -Pop, names_to = "statistic", values_to = "value") 
head(broad_poppr_result)

country_poppr_result = read.csv("populations/country_poppr_results.csv") %>%
  dplyr::select(-c(File, Hexp, SE, E.5, eMLG)) %>%
  dplyr::rename("PA" = "private_alleles",
                "AR" = "allelerichness") %>%
  tidyr::pivot_longer(cols = -Pop, names_to = "statistic", values_to = "value")
head(country_poppr_result)

site_poppr_result = read.csv("populations/site_poppr_results.csv") %>%
  dplyr::select(-c(File, Hexp, SE, E.5, eMLG)) %>%
  dplyr::rename("PA" = "private_alleles",
                "AR" = "allelerichness") %>%
  tidyr::pivot_longer(cols = -Pop, names_to = "statistic", values_to = "value")
head(site_poppr_result)


popstats.broad.plot = ggplot(broad_poppr_result, aes(x = Pop, y = value, fill = Pop)) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.6, , colour = "black") +
  facet_wrap(~statistic, scales = "free_y") +
  scale_fill_manual(values = c("goldenrod", "darkred", "forestgreen")) +  # Custom fill colors
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Population", y = "Value", title = "") ;popstats.broad.plot

ggsave("figures/popstats.broad.png", popstats.broad.plot, width = 6.2, dpi=400,
       height = 5, units = "in")
ggsave("figures/popstats.broad.svg", popstats.broad.plot, width = 6.2, dpi=400,
       height = 5, units = "in")


country_poppr_result$fill_color <- ifelse(grepl("\\bSA\\b", country_poppr_result$Pop), "forestgreen",
                                          ifelse(grepl("\\b(USA|MEX)\\b", country_poppr_result$Pop), "darkred", "goldenrod"))

popstats.country.plot = ggplot(country_poppr_result, aes(x = Pop, y = value, 
                                                         fill = fill_color)) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.6, colour ="black") +
  facet_wrap(~statistic, scales = "free_y", nrow = 4, ncol = 3) +
  scale_fill_identity() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Population", y = "Value", title = "") ;popstats.country.plot

ggsave("figures/popstats.country.png", popstats.country.plot, width = 6.5, dpi=400,
       height = 6.5, units = "in")
ggsave("figures/popstats.country.svg", popstats.country.plot, width = 6.5, dpi=400,
       height = 6.5, units = "in")


site_poppr_result$fill_color = ifelse(grepl("^SA", site_poppr_result$Pop), "forestgreen",
                                      ifelse(grepl("^(USA|MEX)", site_poppr_result$Pop), "darkred", "goldenrod"))

popstats.site.plot = ggplot(site_poppr_result, aes(x = Pop, y = value, 
                                                   fill = fill_color)) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.6, , colour = "black") +
  facet_wrap(~statistic, scales = "free_y", nrow = 4, ncol = 3) +
  scale_fill_identity() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Population", y = "Value", title = "") ;popstats.site.plot

ggsave("figures/popstats.site.png", popstats.site.plot, width = 7.5, dpi=400,
       height = 6.5, units = "in")