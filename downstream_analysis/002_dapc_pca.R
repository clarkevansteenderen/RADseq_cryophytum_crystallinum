library(SNPfiltR)
library(vcfR)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(PopGenReport)
library(ggplot2)

setwd("D:/RADseq/crystallinum")

popmap.broad = read.table("populations/broad_pops.txt", header = TRUE)
popmap.broad

popmap.country = read.table("populations/country_pops.txt", header = TRUE)
popmap.country

popmap.site = read.table("populations/site_pops.txt", header = TRUE)
popmap.site

cryo.snps = vcfR::read.vcfR("populations/populations.snps.filtered_mac3.vcf")

# convert to genind
cryo.snps_genind = vcfR::vcfR2genind(cryo.snps)
cryo.snps_genind

# convert to genlight
cryo.snps_genlight = vcfR::vcfR2genlight(cryo.snps)
cryo.snps_genlight@ind.names

# Extract sample names from the genlight object
filtered_samples = cryo.snps_genlight@ind.names

# Extract original sample names from popmap
original_samples = popmap.broad$id  # assuming popmap.broad$id is a character vector

# Find missing samples
missing_samples = setdiff(original_samples, filtered_samples)

# Print them
cat("Samples present in popmap but missing from genlight:\n")
print(missing_samples)

# Filter out those missing samples from each popmap
# these were the low quality samples filtered out earlier
popmap.broad = popmap.broad[!popmap.broad$id %in% missing_samples, ]
popmap.country = popmap.country[!popmap.country$id %in% missing_samples, ]
popmap.site = popmap.site[!popmap.site$id %in% missing_samples, ]

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

#####################################################################################
# create STRUCTURE files

dartR::gl2faststructure(x = broad_grouping.genlight, outpath=getwd(),
                outfile = "populations/faststructure_input.str")

# create nexus file -> in the file, change integerdata to standard
dartR::gl2structure(x = broad_grouping.genlight, outpath=getwd(),
                    outfile = "populations/structure_input.str")



# save a version for SplitsTree
sample.div = StAMPP::stamppNeisD(broad_grouping.genlight, pop = F)

#export for splitstree
StAMPP::stamppPhylip(distance.mat=sample.div, 
                     file="populations/splitstree_mac3.txt")

###########################################################################
# BROAD GROUPING PCA
##########################################################################

broad_myCol = c("gold", "red", "forestgreen")

###################################
# STANDARD PCA FIRST
###################################

broad.pca = adegenet::glPca(broad_grouping.genlight, nf = 3)

# Extract PCA scores
scores = as.data.frame(broad.pca$scores)  # PCA coordinates (individuals × axes)
scores$group = pop(broad_grouping.genlight)

head(scores)

scores = scores %>%
  dplyr::filter(!rownames(.) %in% c("HKSA1", "HKSA2", "HKSA3"))

pca.ggplot = ggplot(scores, aes(x = PC1, y = PC2, fill = group)) +
  geom_point(alpha = 0.5, size = 1, shape = 21) +
  labs(
    #title = "PCA of RADseq data",
    x = paste0("PC1 (", round(broad.pca$eig[1] / sum(broad.pca$eig) * 100, 1), "%)"),
    y = paste0("PC2 (", round(broad.pca$eig[2] / sum(broad.pca$eig) * 100, 1), "%)")
  ) +
  scale_fill_manual(values = c(
    "introduced" = "gold",
    "invaded" = "red",
    "native" = "forestgreen"
  )) + 
  ggrepel::geom_text_repel(aes(label = rownames(scores)), 
                           size = 2, colour = "grey50",
                           max.overlaps = Inf,
                           segment.size = 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  #ggtitle("d)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.1)

pca.ggplot

ggsave("figures/pca_mac3_SA.svg", pca.ggplot, width = 6, dpi=400,
       height = 4, units = "in")
ggsave("figures/pca_mac3_SA.png", pca.ggplot, width = 6, dpi=500,
       height = 4, units = "in")

###################################
# DAPC
###################################

broad.dapc = adegenet::dapc(broad_grouping.genlight, var.contrib = TRUE, 
                  scale = FALSE, n.pca = 40, n.da = nPop(broad_grouping.genlight) - 1)

ade4::scatter(broad.dapc, clabel = 0, legend = TRUE, cell=0,
              cstar = 0, cex = 4, pch = 20, solid = 0.4, scree.da=FALSE,
              col = broad_myCol, bg="white",posi.pca="bottomleft",
              cleg = 0.75, xax = 1, yax = 2)

# Extract scores and labels
dapc_broad_df = as.data.frame(broad.dapc$ind.coord)
dapc_broad_df$group = popmap.broad$pop
head(dapc_broad_df)

dapc_broad_df$sample = rownames(dapc_broad_df)

# Optional: mean points per group
centroids.broad = aggregate(. ~ group, data = dapc_broad_df, FUN = mean)

# Plot with repelled labels
broad.dapc = ggplot(dapc_broad_df, aes(x = LD1, y = LD2, fill = group)) +
  geom_point(alpha = 0.5, size = 4, shape = 21) +
  #geom_point(data = centroids, aes(x = LD1, y = LD2), 
  #           shape = 21, size = 4, alpha = 0.8) +
  scale_fill_manual(values=broad_myCol) +
  ggrepel::geom_text_repel(aes(label = sample), size = 2, max.overlaps = 100) +
  ggrepel::geom_text_repel(data = centroids.broad, aes(label = group), 
                           size = 3, #fontface = "bold",
                           max.overlaps = 100) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("a)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50");broad.dapc


###########################################################################
# COUNTRY GROUPING PCA
##########################################################################

country_myCol = c("gold", "gold", "gold", "red", "gold", "forestgreen",
                   "red")

country.dapc = adegenet::dapc(country_grouping.genlight, var.contrib = TRUE, 
                              scale = FALSE, n.pca = 40, 
                              n.da = nPop(country_grouping.genlight) - 1)

ade4::scatter(country.dapc, posi.pca = "topleft", clabel = 0, leg = TRUE,
              cstar = 0, cell = 0, solid = 0.6, pch = 20, cex = 3, scree.da=FALSE,
              col = country_myCol)

# Extract scores and labels
dapc_country_df = as.data.frame(country.dapc$ind.coord)
dapc_country_df$group = popmap.country$pop
head(dapc_country_df)

# Optional: mean points per group
centroids.country = aggregate(. ~ group, data = dapc_country_df, FUN = mean)

# Plot with repelled labels
country.dapc = ggplot(dapc_country_df, aes(x = LD1, y = LD2, fill = group)) +
  geom_point(alpha = 0.5, size = 4, shape = 21) +
  #geom_point(data = centroids, aes(x = LD1, y = LD2), 
  #           shape = 21, size = 4, alpha = 0.8) +
  scale_fill_manual(values=country_myCol) +
  ggrepel::geom_text_repel(data = centroids.country, aes(label = group), 
                           size = 3, #fontface = "bold",
                           max.overlaps = 100) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("b)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50");country.dapc

###########################################################################
# SITE GROUPING PCA
##########################################################################

# Initialize with goldenrod
pop_levels = levels(pop(site_grouping.genlight))
site_myCol = rep("goldenrod", length(pop_levels))
names(site_myCol) = pop_levels
# Set acronyms (2-character entries) to forestgreen
site_myCol[nchar(pop_levels) == 2] = "forestgreen"
# Override specific countries with red
site_myCol[c("UnitedStates", "Mexico")] = "red"


nInd(site_grouping.genlight)
length(pop(site_grouping.genlight))
levels(pop(site_grouping.genlight))

# this gave an error when there was just one representative sample for a site,
# here, the Canary Islands. After removing that sample (CI11.1.1), it worked

site.dapc = adegenet::dapc(site_grouping.genlight, var.contrib = TRUE, 
                              scale = FALSE, n.pca = 40, 
                              n.da = nPop(site_grouping.genlight) - 1)

ade4::scatter(site.dapc, posi.pca = "topleft", clabel = 0, leg = TRUE,
              cstar = 0, cell = 0, solid = 0.6, pch = 20, cex = 3, scree.da=FALSE,
              col = site_myCol)

# Extract scores and labels
dapc_site_df = as.data.frame(site.dapc$ind.coord)
dapc_site_df$group = popmap.site$pop
head(dapc_site_df)

# Optional: mean points per group
centroids.site = aggregate(. ~ group, data = dapc_site_df, FUN = mean)

# Plot with repelled labels
site.dapc = ggplot(dapc_site_df, aes(x = LD1, y = LD2, fill = group)) +
  geom_point(alpha = 0.5, size = 4, shape = 21) +
  #geom_point(data = centroids, aes(x = LD1, y = LD2), 
  #           shape = 21, size = 4, alpha = 0.8) +
  scale_fill_manual(values=site_myCol) +
  ggrepel::geom_text_repel(data = centroids.site, aes(label = group), 
                           size = 3, #fontface = "bold",
                           max.overlaps = 100) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("c)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50");site.dapc

gridExtra::grid.arrange(broad.dapc, country.dapc, site.dapc,
                        nrow = 2, ncol = 2)

clustering.plots = gridExtra::grid.arrange(broad.dapc, country.dapc, site.dapc, pca.ggplot,
                        nrow = 2, ncol = 2)

ggsave("figures/clustering.png", clustering.plots, width = 7, dpi=400,
       height = 7, units = "in")
ggsave("figures/clustering.svg", clustering.plots, width = 7, dpi=400,
       height = 7, units = "in")

#################################################################
# MAP OF SA
#################################################################

gps.points = read.csv("sample_info/all_samples.csv")
gps.points$lat = as.numeric(gps.points$lat)
gps.points$lon = as.numeric(gps.points$lon)

southafrica_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name %in% c("South Africa", "Lesotho", "eSwatini"))

# Plot GPS points on world map to check our locality data is correct 
distr_map = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = southafrica_ext, alpha = 0.5, fill = "grey80") +
  # Add GPS points 
  
  ggrepel::geom_text_repel(
    data = gps.points,
    aes(x = lon, y = lat, label = sample_id),  # replace `id` with your sample ID column
    size = 2.5,
    color = "grey15",
    max.overlaps = 50,   # adjust as needed
    box.padding = 0.25,
    point.padding = 0.1,
    segment.color = "grey80"
  ) +
  
  geom_point(
    data = gps.points, 
    aes(x = lon, y = lat),
    size = 4,
    color = "forestgreen",
    alpha = 0.75
  ) +
  
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    xlim = c(16, 25),      
    ylim = c(-35, -28),
    expand = FALSE
  ) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  #ggtitle(expression(italic("Orachrysops niobe") * ", South Africa")) +
  ggspatial::annotation_scale(
    location = "bl",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  theme_classic() +
  ggtitle("e)")

distr_map

ggsave("figures/distr_map.png", distr_map, width = 7, dpi=400,
       height = 7, units = "in")
ggsave("figures/distr_map.svg", distr_map, width = 7, dpi=400,
       height = 7, units = "in")

# Plot GPS points on world map to check our locality data is correct 
SA_map = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = southafrica_ext, alpha = 0.5, fill = "grey80") +
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    xlim = c(16, 33),      
    ylim = c(-35, -22),
    expand = FALSE
  ) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  #ggtitle(expression(italic("Orachrysops niobe") * ", South Africa")) +
  ggspatial::annotation_scale(
    location = "br",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  theme_classic() 

SA_map

ggsave("figures/SA_map.png", SA_map, width = 7, dpi=400,
       height = 7, units = "in")
ggsave("figures/SA_map.svg", SA_map, width = 7, dpi=400,
       height = 7, units = "in")
