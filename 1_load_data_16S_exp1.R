library(phyloseq)
library(readr)
library(dplyr)
library(phytools)

load("/Users/carameyer/Desktop/PhD/exp1repo/global_env.RData")
# Load raw biom file and clean data

# import biom file
tmp_ps = import_biom("Data/otu_16S_exp1.biom")
# remove extra taxa rank
tmp_ps@tax_table <- tmp_ps@tax_table[,1:7]
# set taxa rank names
colnames(tmp_ps@tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# import mapping file
tmp_mapping_file <- read_csv("Data/mapping_file_16S_exp1.csv")
tmp_mapping_file<- tmp_mapping_file[1:200,]
tmp_design = sample_data(tmp_mapping_file)
sample_names(tmp_design) <- tmp_design$id
# import tree
tmp_tree = read.newick("Data/phylogeny_16S.tre")
#merge onto one phyloseq object
ps_16S_raw = merge_phyloseq(tmp_ps,tmp_design,tmp_tree)

# final ps_16S_raw: 200 samples & 4507 OTUs

rm(list = names(.GlobalEnv)[grep("tmp_",names(.GlobalEnv))])


ps_16S = ps_16S_raw
# remove aberrant samples
ps_16S <- prune_taxa(taxa_sums(ps_16S) > 0, ps_16S) # 3186
# remove samples with seq count < 11000
ps_16S <- prune_samples(sample_sums(ps_16S) > 11000, ps_16S ) # 198 samples 3186 OTUs
ps_16S <- prune_taxa(taxa_sums(ps_16S) > 0, ps_16S) # 3184 OTUs


## Filter most abundant OTUs
tmp_ps = ps_16S
# calculate OTU relative abundance
tmp_df_otu <- as.data.frame(otu_table(tmp_ps))
tmp_df_otu_freq <- apply(tmp_df_otu, 2, FUN=function(x) x/sum(x)*100)
# keep only OTUs with a relative abundance of at least 0.5%
tmp <- apply(tmp_df_otu_freq, 1, FUN=function(x) sum(x>(0.5)))
# select OTUs above frequency threshold
tmp_otus_F1 <- rownames(tmp_df_otu[-which(tmp==0),])
# subset selected OTUs
tmp_ps_filter1 <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F1, tmp_ps) # 244 OTUs


# /!\ takes from 30min to 3h /!\

tmp_ps = tmp_ps_filter1 
# calculate OTUs prevalence in treatment (ttt)
tmp_df <- psmelt(tmp_ps)
tmp_otu_prev_ttt <- data.frame(matrix(ncol=length(unique(tmp_df$treatment)),
                                      nrow=length(unique(tmp_df$OTU)), 
                                      dimnames=list(unique(tmp_df$OTU),
                                                    unique(tmp_df$treatment))), check.names = FALSE)
for (i in unique(tmp_df$OTU)) {
  for (j in unique(tmp_df$treatment)) {
    tmp_otu_prev_ttt[i,j] <- sum(tmp_df$Abundance[tmp_df$OTU == i & tmp_df$treatment == j] > 0,
                                 na.rm = T) / length(tmp_df$Sample[tmp_df$OTU == i & tmp_df$treatment == j]) *100
  }
  
} 
rm(i,j)
# calculate maximum OTUs prevalence in treatment
tmp <- apply(tmp_otu_prev_ttt,1, FUN=function(x) max(x))
# select OTUs above a minimum prevalence in treatment threshold set to 60% 
tmp_otus_F2 <- rownames(tmp_otu_prev_ttt[which(tmp >= 60),])
# subset selected OTUs
tmp_ps <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F2, tmp_ps)
ps_16S_most_abund = tmp_ps

# 226 OTUs in 198 samples

rm(list = names(.GlobalEnv)[grep("tmp",names(.GlobalEnv))])

# Try less stringent rel abund filter
tmp_ps = ps_16S
# calculate OTU relative abundance
tmp_df_otu <- as.data.frame(otu_table(tmp_ps))
tmp_df_otu_freq <- apply(tmp_df_otu, 2, FUN=function(x) x/sum(x)*100)
# keep only OTUs with a relative abundance of at least 0.02%
tmp <- apply(tmp_df_otu_freq, 1, FUN=function(x) sum(x>(0.02)))
# select OTUs above frequency threshold
tmp_otus_F1_new <- rownames(tmp_df_otu[-which(tmp==0),])
# subset selected OTUs
tmp_ps_filter1_new <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F1_new, tmp_ps) # 1392 OTUs

# keep only OTUs with a relative abundance of at least 0.2%
tmp <- apply(tmp_df_otu_freq, 1, FUN=function(x) sum(x>(0.2)))
# select OTUs above frequency threshold
tmp_otus_F1.2 <- rownames(tmp_df_otu[-which(tmp==0),])
# subset selected OTUs
tmp_ps_filter1.2 <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F1.2, tmp_ps) # 444 OTUs

# keep only OTUs with a relative abundance of at least 0.1%
tmp <- apply(tmp_df_otu_freq, 1, FUN=function(x) sum(x>(0.1)))
# select OTUs above frequency threshold
tmp_otus_F1.1 <- rownames(tmp_df_otu[-which(tmp==0),])
# subset selected OTUs
tmp_ps_filter1.1 <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F1.1, tmp_ps) # 661 OTUs

# Now filter out based on prevelance in treatments
# /!\ takes from 30min to 3h /!\

tmp_ps = tmp_ps_filter1.1 
# calculate OTUs prevalence in treatment (ttt)
tmp_df <- psmelt(tmp_ps)
tmp_otu_prev_ttt <- data.frame(matrix(ncol=length(unique(tmp_df$treatment)),
                                      nrow=length(unique(tmp_df$OTU)), 
                                      dimnames=list(unique(tmp_df$OTU),
                                                    unique(tmp_df$treatment))), check.names = FALSE)
for (i in unique(tmp_df$OTU)) {
  for (j in unique(tmp_df$treatment)) {
    tmp_otu_prev_ttt[i,j] <- sum(tmp_df$Abundance[tmp_df$OTU == i & tmp_df$treatment == j] > 0,
                                 na.rm = T) / length(tmp_df$Sample[tmp_df$OTU == i & tmp_df$treatment == j]) *100
  }
  
} 
rm(i,j)
# calculate maximum OTUs prevalence in treatment
tmp <- apply(tmp_otu_prev_ttt,1, FUN=function(x) max(x))
# select OTUs above a minimum prevalence in treatment threshold set to 60% 
tmp_otus_F2.1 <- rownames(tmp_otu_prev_ttt[which(tmp >= 60),])
# subset selected OTUs
tmp_ps <- prune_taxa(taxa_names(tmp_ps) %in% tmp_otus_F2.1, tmp_ps)
ps_16S_most_abund.1 = tmp_ps # 568 OTUs in 198 samples


otu_abund_filtered <- as.data.frame(otu_table(tmp_ps_filter1.1))
otu_abund_filtered$total_count <- rowSums(otu_abund_filtered)
sorted <- otu_abund_filtered[order(otu_abund_filtered$total_count),]

least <- row.names(sorted[1:10,])
least

least_rel <- subset(sign.Control.rel, sign.Control.rel$OTU %in% 
                     least) 

prev <- subset(sorted_otus, rownames(sorted_otus) %in% rownames(prev.control_U))
prev <- prev[, c(199,1:198)]

otu_sorted_not_overlapping_U <- subset(sorted_otus, rownames(sorted_otus) %in% not_overlapping_U)
otu_sorted_not_overlapping_U <- otu_sorted_not_overlapping_U[, c(199, 1:10, 50:59, 100:109, 150:159)]

tmp_ps_F1 = subset_samples(ps_16S_most_abund.1,ps_16S_most_abund.1@sam_data$fraction  == "F1")
otu_abund_F1 <- as.data.frame(otu_table(tmp_ps_F1))
otu_abund_F1$total_count <- rowSums(otu_abund_F1)
sorted_otus_F1 <- otu_abund_F1[order(otu_abund_F1$total_count),]
otu_sorted_not_overlapping_F1_any <- subset(sorted_otus_F1, rownames(sorted_otus_F1) %in% not_overlapping_F1_any)


prev_plus.001 <- subset(sorted_otus, rownames(sorted_otus) %in% otu.pres.U)
                        
otu_abund_f2 <- as.data.frame(otu_table(ps_16S_most_abund.1))
otu_abund_f2$total_count <- rowSums(otu_abund_f2)
sorted_otus <- otu_abund_f2[order(otu_abund_f2$total_count),]

least_abund <- row.names(sorted_otus[1:10,])
least_abund
least_abund_f2 <- subset(sign.Control.rel, sign.Control.rel$OTU %in% 
                      least_abund) 
least_abund_f2$min <- apply(least_abund_f2, 1, FUN = function(x) {min(x[x > 1e-06])})
least_abund_f2 <- least_abund_f2[order(least_abund),]
least_abund_f2 <- least_abund_f2[match(least_abund, row.names(least_abund_f2)),]
order(rownames(least_abund_f2)) <- least_abund


ps = prune_taxa(taxa_names(ps_16S_most_abund.1), ps_rel_abund)

otu_rel_f2 <- as.data.frame(otu_table(ps))
otu_rel_f2$sum_rel_abund <- rowSums(otu_rel_f2)
sorted_otus_rel <- otu_rel_f2[order(otu_rel_f2$sum_rel_abund),]
sorted_otus_rel <- otu_rel_f2[match(rownames(sorted_otus), rownames(otu_rel_f2)),]

