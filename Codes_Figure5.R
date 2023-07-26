## A biobank of pediatric patient-derived-xenograft models in cancer precision medicine trial MAPPYACTS for relapsed and refractory tumors
## Figure 5: HLA class I alleles detected in PDX models from pediatric solid cancers
pacman::p_load(readxl, tidyverse, pheatmap, RColorBrewer, ggplot2)

df <- read_excel("Supplementary Data 3.xlsx", sheet = "Figure5") %>% 
  filter(Gene == "HLA-A") %>% # select HLA gene (HLA-A, HLA-B or HLA-C)
  select(-c(Gene)) %>% 
  column_to_rownames("Allele")

## Plot frequencies for controls
controls <- df %>% select(`Frequency in control`, Supertype) %>% rownames_to_column("Allele") %>% 
  mutate(Cohort = "EUR controls")
controls$Supertype <- factor(controls$Supertype, levels = c("A2", "A3", "A1", "A24", "A26", "NA"))
controls$Supertype <- factor(controls$Supertype, levels = c("B44", "B7", "B8", "B39", "B62", "B58", "B27", "NA"))
controls$Supertype <- factor(controls$Supertype, levels = c("C1", "C2"))

p <- ggplot(controls, aes(x = Allele, y = Cohort, size = `Frequency in control`)) +
  geom_point() +
  scale_size_continuous(name = "Allele frequency", range = c(0, 3), limits = c(0, 0.3), 
                        breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25)) +
  scale_y_discrete(position = "right") + 
  theme(axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.line = element_line(), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(), 
        legend.position = "none",
        panel.spacing = unit(0.2, "lines"), 
        strip.text.x = element_text(size = 7),
        strip.background = element_rect(fill = "white", colour = "black")
  ) + facet_grid(~ Supertype, scales = "free", space = "free")

## Plot frequencies for cohorts
df.plot <- df %>% select(-c(`Frequency in control`, Supertype))

# Create a matrix for the numbers to be plotted in the heatmap
df.tmp <- df.plot
df.tmp[df.tmp == 0] <- NA
df.tmp[is.na(df.tmp)] <- ""

breaksList <- seq(0, 12, by = 0.05)
pheatmap(t(df.plot), display_numbers = t(df.tmp), number_format = '%.0f',
         labels_row = c("  OS", "  EWS", "  RMS", "    eRMS", "    aRMS",
                        "  NRSTS", "  NB", "    MYCN-NA", "    MYCN-A", "    11qWT", "    11qLOH",
                        "  NPB", "  HPB", "  HGG", "  EP"),
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(6, 13), gaps_col = c(4, 13, 18, 20, 25),
         angle_col = "90", fontsize = 7, fontsize_number = 8, 
         cellwidth = 10, cellheight = 10,
         color = colorRampPalette(c("white", "#FFFFCC", "#FED976", "#FEB24C",
                                    "#FD8D3C", "#FC4E2A", "#D93E1E"))(length(breaksList)),
         breaks = breaksList
)
## Control- and cohort-specific plots were combined manually by HLA gene
