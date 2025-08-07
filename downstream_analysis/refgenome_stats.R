library(ggplot2)

setwd("D:/RADseq/crystallinum")

# How well did each sample fragments align to the ref genome?

refgen.alignment.stats = read.delim("logfiles/refgenome_alignment_stats.txt",
                                    header = FALSE)

names(refgen.alignment.stats) = c("id", "percentage")
head(refgen.alignment.stats)

summary(refgen.alignment.stats)
Rmisc::summarySE(refgen.alignment.stats, measurevar = "percentage")

refgenstats.plot = ggplot(refgen.alignment.stats, aes(x = id, y = percentage, group = 1)) +
  geom_line(lwd = 0.65, colour = "lightblue") +     # Line plot
  geom_point(shape = 21, size = 1, fill = "darkblue", aes(alpha = 0.2)) +    # Points at each K
  labs(x = "Sample", y = "Percentage") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6)) +
  ggtitle("Reference genome alignment success")

refgenstats.plot

ggsave("figures/refgenome_stats.png", refgenstats.plot, 
       width = 8, height = 4, units = "in")

ggplot(refgen.alignment.stats, aes(x = percentage)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(
    title = "Density Plot of Alignment Percentages",
    x = "Alignment Percentage",
    y = "Density"
  ) + 
  theme_classic() +
  geom_vline(aes(xintercept = mean(percentage)), color = "red", linetype = "dashed")

ggplot(refgen.alignment.stats, aes(x = percentage)) +
  geom_point(
    aes(y = 0),  # Puts all points in a line
    position = position_jitter(height = 0.1),  # Optional: adds slight vertical jitter
    size = 3,
    color = "steelblue"
  ) +
  labs(
    title = "Alignment Percentages as Points",
    x = "Alignment Percentage",
    y = NULL
  ) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


####################################################################
# if we only have the stderr output to work with
####################################################################

log_lines = readLines(
  "populations_guerich_included/refgen_align_stderr.txt")

# Extract lines containing "overall alignment rate"
alignment_lines <- grep("overall alignment rate", log_lines, value = TRUE)

# Extract just the percentage value using a regular expression
alignment_percentages <- sub("^([0-9.]+)% overall alignment rate$", "\\1", alignment_lines)

# Convert to numeric if you want to do further analysis
alignment_percentages <- as.numeric(alignment_percentages)

# Create a data frame
alignment_df <- data.frame(alignment = alignment_percentages)

# now read in the original barcodes file for samples on both plates

bothplates.pops = read.csv("populations_guerich_included/bothplates_pops.txt",
                           header = FALSE, sep = "\t") 
colnames(bothplates.pops) = c("sample_id", "pop")

alignment_df$sample_id = bothplates.pops$sample_id


refgenstats.plot = ggplot(alignment_df, aes(x = sample_id, y = alignment, group = 1)) +
  geom_line(lwd = 0.65, colour = "lightblue") +     # Line plot
  geom_point(shape = 21, size = 1, fill = "darkblue", aes(alpha = 0.2)) +    # Points at each K
  labs(x = "Sample", y = "Percentage") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 5)) +
  ggtitle("Reference genome alignment success")

refgenstats.plot

ggsave("populations_guerich_included/refgenome_stats.png", refgenstats.plot, 
       width = 9, height = 4, units = "in")
