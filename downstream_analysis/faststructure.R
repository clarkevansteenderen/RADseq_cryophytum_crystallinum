library(ggplot2)
library(tidyr)
suppressMessages(library(dplyr))
library(RColorBrewer)
library(magrittr)

setwd("D:/RADseq/crystallinum")
source("R_scripts/downstream_analysis/structure.plot.R")

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# Set the theme for the maps
theme_opts <- list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(colour = "black"),
    panel.border = element_rect(fill = NA),
    legend.key = element_blank()
  )
)

# faststructure output:

#Model complexity that maximizes marginal likelihood = 12
#Model components used to explain structure in data = 4

K2struc.file = "populations/faststructure_output/fastStructure.2.meanQ"
K3struc.file = "populations/faststructure_output/fastStructure.3.meanQ"
K4struc.file = "populations/faststructure_output/fastStructure.4.meanQ"
K5struc.file = "populations/faststructure_output/fastStructure.5.meanQ"
K12struc.file = "populations/faststructure_output/fastStructure.12.meanQ"

sample.info = read.delim("populations/broad_pops_filtered.txt", 
                         header = T, sep = " ")
names(sample.info) = c("id", "pop")

sample.info = sample.info %>%
  mutate(
    pop = case_when(
      pop == "introduced" ~ "Introduced",
      pop == "invaded" ~ "Invaded",
      pop == "native" ~ "Native",
      # Add more countries as needed
      TRUE ~ pop  # Keep the original name if it doesn't match any condition
    )
  )
head(sample.info)

K2 = structure.plot(K2struc.file, sample.info, kval = "2"); K2
K3 = structure.plot(K3struc.file, sample.info, kval = "3"); K3
K4 = structure.plot(K4struc.file, sample.info, kval = "4"); K4
K12 = structure.plot(K12struc.file, sample.info, kval = "12"); K12

strucplots = gridExtra::grid.arrange(K2, K4, ncol=1)
ggsave("figures/structureplot.png", strucplots, width = 10, height = 6, units = "in")
ggsave("figures/structureplot.svg", strucplots, width = 10, height = 6, units = "in")

###################
# MAPS
###################

sample.info.site = read.delim("populations/site_pops_filtered.txt", 
                         header = T, sep = " ")
names(sample.info.site) = c("id", "pop")
head(sample.info.site)

gps.data = read.csv("sample_info/all_samples.csv") %>%
  janitor::clean_names() %>%
  mutate(latitude = ifelse(lat == "Not given", NA, lat),
         longitude = ifelse(lon == "None given", NA, lon)) %>%
  dplyr::select(sample_id, latitude, longitude, origin) %>%
  rename(name = sample_id,
         lat = latitude, lon = longitude)

gps.data$lat = as.numeric(gps.data$lat)
gps.data$lon = as.numeric(gps.data$lon)
gps.data$origin = as.factor(gps.data$origin)
head(gps.data)

# Keep only GPS points that match the sample IDs in sample.info.site
gps.data.pruned <- gps.data %>%
  filter(name %in% sample.info.site$id)

K4hapmap = hapmap.global(K4struc.file, sample.info.site, gps.data.pruned,
                         piesize = 0.06) +
  coord_sf(
    crs = 4326,
    xlim = c(16, 25),      
    ylim = c(-35, -28),
    expand = FALSE
  ) +
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
  );K4hapmap

ggsave("figures/hapmap.png", K4hapmap, width = 10, height = 6, units = "in", dpi = 400)
ggsave("figures/hapmap.svg", K4hapmap, width = 10, height = 6, units = "in", dpi = 400)


#######################################
# global distribution
#######################################

world_map = rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

global_distr = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = world_map, alpha = 0.5) +
  # Add GPS points 
  geom_point(
    data = na.omit(gps.data.pruned), 
    size = 2, 
    aes(x = lon, y = lat
        #fill = sample_status, shape = sample_status
        ),  # Aesthetic for fill and shape
    #color = "black",  # Outline color of points
    stroke = 0.6  # Stroke to ensure the shape outlines are visible
  ) +
  # Set custom colors for the fill based on sample_status
  #scale_fill_manual(values = c("forestgreen", "darkorange", "darkred")) +  
  # Set custom shapes for each sample_status
  #scale_shape_manual(values = c(21, 24, 22)) +  
  coord_sf(
    crs = "+proj=eqearth",
    expand = FALSE,
    ylim = c(-40, 50),
    xlim = c(-130, 50)
  ) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.03, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  theme(legend.position = "top") +
  labs(fill = "Status")   # Label for the color legend
  #guides(
    #fill = guide_legend(override.aes = list(shape = c(21, 24, 22), 
                                        #    fill = c("forestgreen", "darkorange", "darkred"))),  # Combine fill and shape
    #shape = "none"  # Remove the shape legend (since it's included in the fill legend)
  #)

global_distr

ggsave("figures/world_distrib.png", global_distr, height=4, width=6)
ggsave("figures/world_distrib.svg", global_distr, height=4, width=6)
