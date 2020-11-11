rm(list = ls())
setwd("~/git_repos/mink_homoplasies")
require(tidyverse); require(ggplot2); require(data.table); require(reshape2); require(ggpubr)
mink <- read.csv("data/Combined_05_11_2020.deMaio.BASECOUNTS.parsed.csv", 
                       stringsAsFactors = F,
                       check.names = F)

human <- read.csv("../covid_homoplasy_2/gisaid_cov2020_sequences.30.07.2020.QC.human.nextstrain_filter.QC.deMaiomask.HPfinder.noambig.BASECOUNTS.csv",
                  stringsAsFactors = F,
                  check.names = F)

# parse <- function(df) {
#   colnames(df)[1] <- "ref"
#   df <- df[, c("ref", "a", "c", "g", "t")]
#   colnames(df) <- c("ref", "A", "C", "G", "U")
#   df$ref <- toupper(df$ref)
#   df$ref[df$ref == "T"] <- "U"
#   mut_list <- list()
#   
#   for (ref in c("A", "C", "U", "G")) {
#     for (variant in c("A", "C", "U", "G")){
#       if (ref != variant) {
#         morsel <- df[df$ref == ref, variant]
#         mut_list[[paste0(ref, "->", variant)]] <- sum(morsel)
#       }
#     }
#   }
#   
#   intermediate_df <- reshape2::melt(mut_list)
#   mut_types <- c("C->U", "A->G", "G->U", "C->A", "A->U", "G->C")
#   mut_rev <- c("U->C", "G->A", "U->G", "A->C", "U->A", "C->G")
#   plot_df <- data.frame(mutation = mut_types, value = 0)
#   
#   for (i in seq(length(mut_types))) {
#     mut <- mut_types[i]
#     rev <- mut_rev[i]
#     ratio <- intermediate_df[intermediate_df$L1 == mut, ]$value / 
#       intermediate_df[intermediate_df$L1 == rev, ]$value
#     plot_df[plot_df$mutation == mut, "value"] <- ratio
#   }
#   plot_df$mutation <- as.character(plot_df$mutation)
#   plot_df$mutation <- factor(plot_df$mutation, levels = plot_df$mutation[order(plot_df$value, decreasing = T)])
#   
#   return(list(plot_df, intermediate_df))
# }
# 
# mink_df <- parse(mink)[[1]]
# human_df <- parse(human)[[1]]
# 
# parse(mink)[[2]]
# parse(human)[[2]]
# # Set order of bars
# human_df$mutation <- factor(human_df$mutation, levels = human_df$mutation[order(human_df$value, decreasing = T)])
# mink_df$mutation <- factor(mink_df$mutation, levels = levels(human_df$mutation))
# 
# fills <- c("tomato", "steelblue", "darkseagreen", 
#            "grey", "orange", "firebrick")
# # Plot
# plot_figure <- function(df) {
#   plt <- ggplot(df, aes(x = mutation, y = value, fill = mutation)) +
#     geom_bar(stat = "identity") +
#     geom_text(aes(label = round(value, 2)), position = position_dodge(width=0.9), vjust=0, size = 3) +
#     labs(x = "Mutation Type", y = "Ratio") +
#     theme(axis.text.x = element_text(angle = 30),
#           legend.position = "None") +
#     scale_fill_manual(values = fills)
#   return(plt)
# }
# 
# ggarrange(plot_figure(human_df), plot_figure(mink_df))
# ggsave("snp_counts/mutation_asymmetry.png", dpi = 300, height = 4, width = 8)
############ Plot mutation frequencies ##########################################

parse_indiv <- function(df) {
  colnames(df)[1] <- "ref"
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
  
  plot_df <- reshape2::melt(mut_list)
  colnames(plot_df)[2] <- "mutation"
  plot_df$mutation <- as.character(plot_df$mutation)
  # Proportions
  # plot_df$value <- plot_df$value / sum(plot_df$value)
  
  return(plot_df)
}

mink_df <- parse_indiv(mink)
human_df <- parse_indiv(human)

# Set order of bars
human_df$mutation <- factor(human_df$mutation, levels = human_df$mutation[order(human_df$value, decreasing = T)])
mink_df$mutation <- factor(mink_df$mutation, levels = levels(human_df$mutation))

fills <- c("tomato", "steelblue", "darkseagreen", 
           "grey", "orange", "firebrick", "slategrey", "tomato3", "cyan3", "burlywood",
           "darkgoldenrod3", "darkturquoise")
# Plot
plot_figure <- function(df) {
  plt <- ggplot(df, aes(x = mutation, y = value, fill = mutation)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(width=0.9), vjust=0, size = 2) +
    labs(x = "Mutation Type", y = "Ratio") +
    theme(axis.text.x = element_text(angle = 30),
          axis.title = element_blank(),
          legend.position = "None") +
    scale_fill_manual(values = fills)
  return(plt)
}

# Perform Chi square test
merged_df <- inner_join(mink_df, human_df, by = "mutation", suffix = c("_mink", "_human"))
rownames(merged_df) <- merged_df$mutation
merged_df <- merged_df[, -2]
test <- chisq.test(merged_df, simulate.p.value = T)
p <- round(as.numeric(test$p.value), 6)
chisq <- as.numeric(test$statistic)
  
plt <- ggarrange(plot_figure(human_df), plot_figure(mink_df))
annotate_figure(plt, bottom = "Mutation Type", left = "Frequency")
ggsave("snp_counts/human_mink_snp_counts.png", dpi = 300, height = 4, width = 8)
