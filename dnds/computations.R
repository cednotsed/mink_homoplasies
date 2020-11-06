setwd("~/git_repos/covid_homoplasy_2/")
df <- read.csv("dnds/parsed_savs_n46721_040820.csv", check.names = F, stringsAsFactors = F)

# total no. of C->U
sum(df[df$ref == 'C' & df$snp == 'T', ]$snp_count)
# Total no. of SNPs
sum(df$snp_count)
syn <- df[df$NS_S == "S", "snp_count"]
ns <- df[df$NS_S == "NS", "snp_count"]
# NS/S
sum(ns)/sum(syn)
# NS/(S + NS)
sum(ns) / (sum(syn) + sum(ns))
ns <- df[df$NS_S == "NS", ]
# (NS and C->U) / NS 
sum(ns[ns$ref == 'C' & ns$snp == 'T', ]$snp_count) / sum(ns$snp_count)
CU <- df[df$ref == 'C' & df$snp == 'T', ]
# (NS and C->U) / C->U
sum(CU[CU$NS_S == "NS", "snp_count"]) / sum(CU$snp_count)

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
