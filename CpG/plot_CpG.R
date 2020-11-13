rm(list = ls())
setwd("~/git_repos/mink_homoplasies")
require(tidyverse); require(data.table); require(ggplot2); require(Biostrings); require(see); require(ggpubr)
# fasta <- read.fasta(file = "CpG/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln")
# fasta.subset <- fasta[1:10]
# write.fasta(fasta.subset, names(fasta.subset), file.out = "CpG/human_subset.fasta")

parse <- function(file_path, prefix) {
  fasta <- readDNAStringSet(file = file_path)
  dn <- dinucleotideFrequency(fasta, step=1,
                        as.prob=F, as.matrix=FALSE,
                        fast.moving.side="right", with.labels=TRUE)
  cg <- dn[, "CG"]
  
  return(tibble(cg_freq = cg, host = prefix))
}

mink_cg <- parse("CpG/Combined_05_11_2020.deMaio.rename.aln", "mink")
human_cg <- parse("CpG/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln", "human")

plot_df <- rbind(mink_cg, human_cg) %>%
  mutate(host = recode(host, "mink" = "Mink", "human" = "Human"))

# Get statistics
stat_df <- plot_df %>%
  group_by(host) %>%
  summarise(mean = mean(cg_freq), median = median(cg_freq))

plt <- ggplot(plot_df, aes(x = host, y = cg_freq, fill = host)) +
       geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
       geom_point(position = position_jitter(width = 0.08), 
                  size = 0.5, 
                  alpha = 0.3,
                  color = "black") +
       geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
                    width = 0.2, 
                    outlier.shape = NA, 
                    alpha = 1) +
       labs(y = "CpG Frequency", x = "Host") +
       coord_flip() +
       theme(legend.position = "none",
             plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
       stat_compare_means(method = "wilcox.test", 
                          label.x = 0.5, label.y = 363)

plt
ggsave("CpG/CpG_frequencies.png", plot = plt, dpi = 300, width = 8, height = 5)
