rm(list = ls())
setwd("~/git_repos/covid_homoplasy_2")
df <- read.csv("gisaid_cov2020_sequences.30.07.2020.QC.human.nextstrain_filter.QC.deMaiomask.HPfinder.noambig.BASECOUNTS.csv", 
               stringsAsFactors = F,
               check.names = F)
colnames(df)[1] <- "ref"
# df <- df[266:29674, ]
# exclude <- c(seq(21556, 21562), seq(25385, 25392), seq(26221, 26244), seq(26473, 26522), seq(27192, 27201), seq(27388, 27393), seq(27888, 27893), seq(28260, 28273), seq(29534, 29557))
# df <- df[-exclude, ]
df <- df[, c("ref", "a", "c", "g", "t")]
colnames(df) <- c("ref", "A", "C", "G", "U")
df$ref <- toupper(df$ref)
df$ref[df$ref == "T"] <- "U"
mut_list <- list()

for (ref in c("A", "C", "U", "G")) {
  for (variant in c("A", "C", "U", "G")){
    if (ref != variant) {
      morsel <- df[df$ref == ref, variant]
      mut_list[[paste0(ref, "->", variant)]] <- sum(morsel)
    }
  }
}

require(reshape2)
plot_df <- melt(mut_list)
plot_df$L1 <- as.character(plot_df$L1)
plot_df$L1 <- factor(plot_df$L1, levels = plot_df$L1[order(plot_df$value, decreasing = T)])

# Perform computations
Cs <- plot_df[grep("C->", plot_df$L1), ]
sum(Cs$value)
Cs$value[Cs$L1 == "C->U"]/sum(Cs$value)

# Proportions
plot_df$value <- plot_df$value / sum(plot_df$value)

require(ggplot2)
plt1 <- ggplot(plot_df, aes(x = L1, y = value, fill = L1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(value, 2)), position = position_dodge(width=0.9), vjust=0, size = 3) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 30),
        legend.position = "None") +
  scale_fill_discrete()

# ggsave("snp_counts/all_SNPs_by_pos.png", dpi = 600)

###########################################################################
# Homoplasy positions
meta <- read.csv("Demaio_filtered-homoplasic-sites-table_alt2_vs_alt1_0.2.csv")
meta <- meta$bp

df2 <- df[meta, ]
mut_list <- list()

for (ref in c("A", "C", "U", "G")) {
  for (variant in c("A", "C", "U", "G")){
    if (ref != variant) {
      morsel <- df2[df2$ref == ref, variant]
      mut_list[[paste0(ref, "->", variant)]] <- sum(morsel)
    }
  }
}

plot_df2 <- melt(mut_list)
plot_df2$L1 <- factor(plot_df2$L1, levels = plot_df2$L1[order(plot_df2$value, decreasing = T)])

# Proportions
plot_df2$value <- plot_df2$value / sum(plot_df2$value)

require(ggplot2)
plt2 <- ggplot(plot_df2, aes(x = L1, y = value, fill = L1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(value, 2)), position = position_dodge(width=0.9), vjust=0, size = 3) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 30),
        legend.position = "None") +
  scale_fill_discrete()


require(ggpubr)
figure <- ggarrange(plt1, plt2, ncol = 2, nrow = 1, labels = "auto")
annotate_figure(figure, 
                bottom = text_grob("Mutation Type"),
                left = text_grob("Proportion", rot = 90)
                )
ggsave("snp_counts/snps_by_pos.png", dpi = 600, width = 8, height = 4)

