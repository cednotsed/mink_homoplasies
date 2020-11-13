setwd("~/git_repos/mink_homoplasies/")
df <- read.csv("dnds/parsed_savs_minks_040820.csv", check.names = F, stringsAsFactors = F)
description <- c()
statistic <- c()

# Total no. of C->T
ct <- sum(df[df$ref == 'C' & df$snp == 'T', ]$snp_count)
description <- c(description, "total no. of C->U")
statistic <- c(statistic, ct)

# Total no. of SNPs
total <- sum(df$snp_count)
description <- c(description, "Total no. of SNPs")
statistic <- c(statistic, total)

# No. of synonymous mutations
syn <- sum(df[df$NS_S == "S", "snp_count"])
description <- c(description, "No. of synonymous mutations")
statistic <- c(statistic, syn)

# No. of non-synonymous mutations
non_syn <- sum(df[df$NS_S == "NS", "snp_count"])
description <- c(description, "No. of non-synonymous mutations")
statistic <- c(statistic, non_syn)

# NS/S
ns_s <- non_syn / syn
description <- c(description, "NS/S")
statistic <- c(statistic, ns_s)

# NS/(S + NS)
ns_s_ns <- non_syn / (syn + non_syn)
description <- c(description, "NS/(S + NS)")
statistic <- c(statistic, ns_s_ns)

# (NS and C->U) / NS 
ns <- df[df$NS_S == "NS", ]
ns_cu_ns <- sum(ns[ns$ref == 'C' & ns$snp == 'T', ]$snp_count) / sum(ns$snp_count)
description <- c(description, "(NS and C->U) / NS ")
statistic <- c(statistic, ns_cu_ns)

# (NS and C->U) / C->U
CU <- df[df$ref == 'C' & df$snp == 'T', ]
ns_cu_cu <- sum(CU[CU$NS_S == "NS", "snp_count"]) / sum(CU$snp_count)
description <- c(description, "(NS and C->U) / C->U")
statistic <- c(statistic, ns_cu_cu)

# NS and C->U / NS and C->?
CU <- df[df$ref == 'C' & df$snp == 'T', ]
CX <- df[df$ref == 'C', ]
ns_cu_cx <- sum(CU[CU$NS_S == "NS", "snp_count"]) / sum(CX[CX$NS_S == "NS", "snp_count"])
description <- c(description, "NS and C->U / NS and C->?")
statistic <- c(statistic, ns_cu_cx)

# all <- read.csv("filtered-homoplasic-sites-table_alt2_vs_alt1_0.2.csv", check.names = F, stringsAsFactors = F)
# all <- all[, c("bp", "N.isolates.with.homoplasy")]
# colnames(all)[2] <- "snp_count"
# merge_df <- merge(all, df)
# dim(merge_df)
# nrow(all[!(all$bp %in% df$bp), ])
# 
# meta <- read.csv("dnds/orf_positions.csv", stringsAsFactors = F)
# meta <- meta[-c(1, 14), ]
# include <- c()
# for (i in seq(nrow(meta))) {
#   include <- c(include, seq(meta[i, "start"], meta[i, "end"]))  
# }
# 
# length(all$bp[!(all$bp %in% include)])
# all[!(all$bp %in% include), ]
# meta

results_df <- data.frame(description = description, statistic = statistic)

fwrite(results_df, "dnds/dnds_results.csv")
