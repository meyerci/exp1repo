library(export)
library(UpSetR)
library(dplyr)
library(tidyr)
library(grid)
library(ggplot2)

# Find significant effects of pesticide in fractions ----
# subset significant changes from differential abundance analysis; add column identifying fraction
sig.otu.U <- subset(glmPLN.pairwise.global.U, p.adjust < 0.05) %>%
  mutate(fraction = "Unfiltered") %>%
  mutate(comparison = paste(fraction, contrast))
sig.otu.F1 <- subset(glmPLN.pairwise.global.F1, p.adjust < 0.05) %>%
  mutate(fraction = "F1") %>%
  mutate(comparison = paste(fraction, contrast))
sig.otu.F2 <- subset(glmPLN.pairwise.global.F2, p.adjust < 0.05) %>%
  mutate(fraction = "F2") %>%
  mutate(comparison = paste(fraction, contrast))
sig.otu.F3 <- subset(glmPLN.pairwise.global.F3, p.adjust < 0.05) %>%
  mutate(fraction = "F3") %>%
  mutate(comparison = paste(fraction, contrast))
sig.otu.F4 <- subset(glmPLN.pairwise.global.F4, p.adjust < 0.05) %>%
  mutate(fraction = "F4") %>%
  mutate(comparison = paste(fraction, contrast))

sig.otu.all <- rbind(sig.otu.U, sig.otu.F1, sig.otu.F2, sig.otu.F3, sig.otu.F4)
sig_otus <- unique(sig.otu.all$OTU) #473

# count the number of OTUs per comparison 
nb.OTU.comp.all <- sig.otu.all %>%
  group_by(comparison) %>%
  summarise(OTU = list(unique(OTU)))
head(nb.OTU.comp.all)
# Significant estimates per OTU in each comparison 
sign.Control.all <- sig.otu.all %>%
  group_by(comparison, OTU) %>%
  summarise(estimates = estimate)
head(sign.Control.all) 

sign.Control.all.df <- tibble(sign.Control.all)

#transform from long to wide tibble
sign.Control.all.df<- sign.Control.all.df %>%
  spread(comparison, estimates)
sign.Control.all.df[is.na(sign.Control.all.df)] <- 0
head(sign.Control.all.df)
dim(sign.Control.all.df)

#write.table(x=sign.Control.E, file="sign.Control.Estimates.csv")

##how many OTUs significant per comparison
nb.OTUs.sign.all <- colSums(sign.Control.all.df != 0) %>% as.data.frame()
nb.OTUs.sign.all
write.table(x=nb.OTUs.sign.all, file="out/nb.OTUs.significant.csv", sep = ";")

###Add % -> % significant OTUs of most abundant OTUs per comparison
nb.OTUs.sign.all.perc <- (colSums(sign.Control.all.df != 0)/ nrow(tax_table(ps_16S_most_abund.1)))*100 
nb.OTUs.sign.all.perc <- as.data.frame(nb.OTUs.sign.all.perc)
nb.OTUs.sign.all.perc
write.table(x=nb.OTUs.sign.all.perc, file="out/nb.OTUs.significant.all.perc.csv", sep = ";")

# Make into a data frame for saving UpSet results later
affected.otus <- data.frame("Affected OTUs" = nb.OTUs.sign.all$., "Percentage of most abundant OTUs" = nb.OTUs.sign.all.perc$nb.OTUs.sign.all.perc,
                            row.names = row.names(nb.OTUs.sign.all.perc))
rownames(affected.otus)[1] <- "Any Treatment"

# what if we make the condition that they have to be present in 60% of the replicates of any treatment in that fraction?
prev.any_F1 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F1 >= 60 | tmp_otu_prev_ttt$`1x_F1` >= 60 |
                            tmp_otu_prev_ttt$`10x_F1` >= 60 | tmp_otu_prev_ttt$`100x_F1` >= 60)
length(unique(rownames(prev.any_F1))) #517
otu.prev.any.F1 <- unique(rownames(prev.any_F1))
length(intersect(unique(sig.otu.F1$OTU), unique(rownames(prev.any_F1)))) #329
length(unique(sig.otu.F1$OTU)) #329 so is it all of them now?
not_overlapping_F1_any <- setdiff(unique(rownames(prev.any_F1)), unique(sig.otu.F1$OTU))

prev.any_F2 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F2 >= 60 | tmp_otu_prev_ttt$`1x_F2` >= 60 |
                        tmp_otu_prev_ttt$`10x_F2` >= 60 | tmp_otu_prev_ttt$`100x_F2` >= 60)
length(unique(rownames(prev.any_F2))) #359

prev.any_F3 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F3 >= 60 | tmp_otu_prev_ttt$`1x_F3` >= 60 |
                        tmp_otu_prev_ttt$`10x_F3` >= 60 | tmp_otu_prev_ttt$`100x_F3` >= 60)
length(unique(rownames(prev.any_F3))) #494

prev.any_F4 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F4 >= 60 | tmp_otu_prev_ttt$`1x_F4` >= 60 |
                        tmp_otu_prev_ttt$`10x_F4` >= 60 | tmp_otu_prev_ttt$`100x_F4` >= 60)
length(unique(rownames(prev.any_F4))) #367

prev.any_U <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_Unfiltered >= 60 | tmp_otu_prev_ttt$`1x_Unfiltered` >= 60 |
                        tmp_otu_prev_ttt$`10x_Unfiltered` >= 60 | tmp_otu_prev_ttt$`100x_Unfiltered` >= 60)
length(unique(rownames(prev.any_U))) #530


all.pres <- left_join("OTU" = all.df$OTU, "present U" = ifelse(otu.prev.any.F1 %in% all.df$OTU, 1, 0))
all.pres <- data.frame(matrix(ncol = 5, nrow = 473))
colnames(all.pres) <- c("U", "F1", "F2", "F3", "F4")
rownames(all.pres) <- all.df$OTU
all.pres[,1] <- apply(all.pres[, 1], 1, function(x) ifelse(rownames(x) %in% rownames(prev.any_U), 1, 0))
all.pres[,1] 
all.pres$U <- ifelse(rownames(all.pres) %in% rownames(prev.any_U), 1, 0)
all.pres$F1 <- ifelse(rownames(all.pres) %in% rownames(prev.any_F1), 1, 0)

prev.control_F2 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F2 >= 60) #from load_data script
length(unique(rownames(prev.control_F2))) #270
otu.pres.F2 <- intersect(unique(rownames(prev.control_F2)), unique(sum.present_F2$OTU))
length(otu.pres.F2) #244

prev.control_F3 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F3 >= 60) #from load_data script
length(unique(rownames(prev.control_F3))) #412
otu.pres.F3 <- intersect(unique(rownames(prev.control_F3)), unique(sum.present_F3$OTU))
length(otu.pres.F3) #378

prev.control_F4 <- subset(tmp_otu_prev_ttt, tmp_otu_prev_ttt$Control_F4 >= 60) #from load_data script
length(unique(rownames(prev.control_F4))) #319
otu.pres.F4 <- intersect(unique(rownames(prev.control_F4)), unique(sum.present_F4$OTU))
length(otu.pres.F4) #294

# Make subsets for OTUs increasing or decreasing in abundance relative to control ----
## Because the contrast is from for example "Control - 1x" a positive estimate means that that OTU is more abundant in the control than 1x
## This means that a negative estimate means increasing in abundance in 1x relative to control
## And positive estimate means decreasing in abundance in 1x relative to control

head(sign.Control.all.df)
#sign.Control.E <- as_data_frame(sign.Control.E)

# make up and down for each DOSE (negative estimates are increasing relative to control)
all.df = as.data.frame(sign.Control.all.df, rownames = 1)
str(all.df)

## UNDER CONSTRUCTION ----
colnames(all.df[2:16])
str(all.df)

all.up

F1_Control_100x <- data.frame(increasing = apply(up.1x[, c(2:6)], 2, function(x) ifelse(x >= 0, 0, 1)))

colnames(all.df) <- c("OTU","F1_control_100x", "F1_control_10x", "F1_control_1x", "F2_control_100x", "F2_control_10x", "F2_control_1x",
                      "F3_control_100x", "F3_control_10x", "F3_control_1x", "F4_control_100x", "F4_control_10x", "F4_control_1x",
                      "Unfiltered_control_100x", "Unfiltered_control_10x", "Unfiltered_control_1x")

#get all OTUs increasing and all decreasing
all.up = all.df
all.up[, -1] <- apply(all.up[, -1], 2, function(x) ifelse(x >= 0, 0, 1))

all.do = all.df
all.do[, -1] <- apply(all.do[, -1], 2, function(x) ifelse(x <= 0, 0, 1))

nb.OTUs.up <- subset(all.up, rowSums(all.ind.up[,-1]) != 0) # 332

nb.OTUs.do <- subset(all.do, rowSums(all.ind.do[,-1]) != 0) # 261

#subset by dose
all.1x <- all.df[ , grep("OTU|1x", colnames(all.df))]
all.10x <- all.df[ , grep("OTU|10x", colnames(all.df))]
all.100x <- all.df[ , grep("OTU|100x", colnames(all.df))]

# make data frame where if OTU is not changing in contrast, it is coded as 1 and if changing is 0
nochange.1x = all.1x
nochange.1x[, -1] <- apply(nochange.1x[, -1], 2, function(x) ifelse(x == 0, 1, 0))

nochange.10x = all.10x
nochange.10x[, -1] <- apply(nochange.10x[, -1], 2, function(x) ifelse(x == 0, 1, 0))

nochange.100x = all.100x
nochange.100x[, -1] <- apply(nochange.100x[, -1], 2, function(x) ifelse(x == 0, 1, 0))

# add condition that OTU must be present in that fraction
# don't have to do anything to pres.nochange.1x$Unfiltered_control_1x because we are assuming all OTUs are present
pres.nochange.1x = nochange.1x
pres.nochange.1x$F1_control_1x <- ifelse(nochange.1x$F1_control_1x == 1 & nochange.1x$OTU %in% rownames(prev.any_F1), 1, 0)
pres.nochange.1x$F2_control_1x <- ifelse(nochange.1x$F2_control_1x == 1 & nochange.1x$OTU %in% rownames(prev.any_F2), 1, 0)
pres.nochange.1x$F3_control_1x <- ifelse(nochange.1x$F3_control_1x == 1 & nochange.1x$OTU %in% rownames(prev.any_F3), 1, 0)
pres.nochange.1x$F4_control_1x <- ifelse(nochange.1x$F4_control_1x == 1 & nochange.1x$OTU %in% rownames(prev.any_F4), 1, 0)

pres.nochange.10x = nochange.10x
pres.nochange.10x$F1_control_10x <- ifelse(nochange.10x$F1_control_10x == 1 & nochange.10x$OTU %in% rownames(prev.any_F1), 1, 0)
pres.nochange.10x$F2_control_10x <- ifelse(nochange.10x$F2_control_10x == 1 & nochange.10x$OTU %in% rownames(prev.any_F2), 1, 0)
pres.nochange.10x$F3_control_10x <- ifelse(nochange.10x$F3_control_10x == 1 & nochange.10x$OTU %in% rownames(prev.any_F3), 1, 0)
pres.nochange.10x$F4_control_10x <- ifelse(nochange.10x$F4_control_10x == 1 & nochange.10x$OTU %in% rownames(prev.any_F4), 1, 0)

pres.nochange.100x = nochange.100x
pres.nochange.100x$F1_control_100x <- ifelse(nochange.100x$F1_control_100x == 1 & nochange.100x$OTU %in% rownames(prev.any_F1), 1, 0)
pres.nochange.100x$F2_control_100x <- ifelse(nochange.100x$F2_control_100x == 1 & nochange.100x$OTU %in% rownames(prev.any_F2), 1, 0)
pres.nochange.100x$F3_control_100x <- ifelse(nochange.100x$F3_control_100x == 1 & nochange.100x$OTU %in% rownames(prev.any_F3), 1, 0)
pres.nochange.100x$F4_control_100x <- ifelse(nochange.100x$F4_control_100x == 1 & nochange.100x$OTU %in% rownames(prev.any_F4), 1, 0)

# For each of 15 UpSets ----
## F1 Control - 100x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F1_100x = all.100x
opp_F1_100x_U <- opp_F1_100x[which(opp_F1_100x$F1_control_100x != 0 & opp_F1_100x$Unfiltered_control_100x != 0),]
opp_F1_100x_U <- opp_F1_100x_U[!(sign(opp_F1_100x_U$F1_control_100x) == sign(opp_F1_100x_U$Unfiltered_control_100x)),]

opp_F1_100x_U <- all.100x %>%
  subset(F1_control_100x != 0 & Unfiltered_control_100x != 0) %>%
  subset(!(sign(F1_control_100x) == sign(Unfiltered_control_100x)))

opp_F1_100x_F2 <- all.100x %>%
  subset(F1_control_100x != 0 & F2_control_100x != 0) %>%
  subset(!(sign(F1_control_100x) == sign(F2_control_100x))) 

opp_F1_100x_F3 <- all.100x %>%
  subset(F1_control_100x != 0 & F3_control_100x != 0) %>%
  subset(!(sign(F1_control_100x) == sign(F3_control_100x)))

opp_F1_100x_F4 <- all.100x %>%
  subset(F1_control_100x != 0 & F4_control_100x != 0) %>%
  subset(!(sign(F1_control_100x) == sign(F4_control_100x)))
opp_F1_100x_F4$OTU

all_opp_F1_100x = all.100x[,-2]
all_opp_F1_100x$Unfiltered_control_100x <- ifelse(all_opp_F1_100x$OTU %in% opp_F1_100x_U$OTU, 1, 0)
all_opp_F1_100x$F2_control_100x <- ifelse(all_opp_F1_100x$OTU %in% opp_F1_100x_F2$OTU, 1, 0)
all_opp_F1_100x$F3_control_100x <- ifelse(all_opp_F1_100x$OTU %in% opp_F1_100x_F3$OTU, 1, 0)
all_opp_F1_100x$F4_control_100x <- ifelse(all_opp_F1_100x$OTU %in% opp_F1_100x_F4$OTU, 1, 0)

diff_F1_100x = pres.nochange.100x
diff_F1_100x$Unfiltered_control_100x <- ifelse(diff_F1_100x$Unfiltered_control_100x == 1 | diff_F1_100x$OTU %in% opp_F1_100x_U$OTU, 1, 0)
diff_F1_100x$F2_control_100x <- ifelse(diff_F1_100x$F2_control_100x == 1 | diff_F1_100x$OTU %in% opp_F1_100x_F2$OTU, 1, 0)
diff_F1_100x$F3_control_100x <- ifelse(diff_F1_100x$F3_control_100x == 1 | diff_F1_100x$OTU %in% opp_F1_100x_F3$OTU, 1, 0)
diff_F1_100x$F4_control_100x <- ifelse(diff_F1_100x$F4_control_100x == 1 | diff_F1_100x$OTU %in% opp_F1_100x_F4$OTU, 1, 0)


F1_control_100x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F1_control_100x, decreasing = all.do$F1_control_100x,
                                    "Unfiltered Different" = diff_F1_100x$Unfiltered_control_100x,
                                    "F2 Different" = diff_F1_100x$F2_control_100x,
                                    "F3 Different" = diff_F1_100x$F3_control_100x,
                                    "F4 Different" = diff_F1_100x$F4_control_100x)
F1_100x_check <- subset(F1_control_100x_check, increasing + decreasing != 0)

F1.100x.upset <- F1_100x_check %>% 
  upset(., sets = c("F4.Different", "F3.Different", "F2.Different", "Unfiltered.Different", "decreasing", "increasing"), keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F1.100x.upset
grid.text(paste(length(F1_100x_check$OTU), "OTUs significantly affected in Fraction 1 100x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_100x_check[which(F1_100x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_100x_check[which(F1_100x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.100x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F1.100x$Increasing.in.all.fractions/affected.F1.100x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.100x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F1.100x$Decreasing.in.all.fractions/affected.F1.100x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F1.100x$Indirect, "of", length(F1_100x_check$OTU), "OTUs (", round(affected.F1.100x$Indirect/length(F1_100x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 


graph2ppt(x= F1.100x.upset ,file="out/F1.100x.upset.pptx")

nrow(F1_100x_check) #255
nrow(F1_100x_check[which(F1_100x_check$increasing == 1),]) #112
nrow(F1_100x_check[which(F1_100x_check$decreasing == 1),]) #143



# count how many direct
otus.F1.100x.up.direct <- F1_100x_check$OTU[which(rowSums(F1_100x_check[,-1]) == 1 & F1_100x_check$increasing == 1)] #44
otus.F1.100x.do.direct <- F1_100x_check$OTU[which(rowSums(F1_100x_check[,-1]) == 1 & F1_100x_check$decreasing == 1)] #22
otus.F1.100x.all.direct <- c(otus.F1.100x.up.direct, otus.F1.100x.do.direct) #66

#vector of indirect otus
`%!in%` <- Negate(`%in%`)
otus.F1.100x.indirect <- F1_100x_check$OTU[which(F1_100x_check$OTU %!in% otus.F1.100x.all.direct)] #189

# count opposite
otus.F1.100x.opp <- unique(c(opp_F1_100x_U$OTU, opp_F1_100x_F2$OTU, opp_F1_100x_F3$OTU, opp_F1_100x_F4$OTU)) #25

# put into data frame to add to affected.otus
affected.F1.100x <- data.frame("Indirect" = length(otus.F1.100x.indirect),
                               "OTUs Increasing" = nrow(F1_100x_check[which(F1_100x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(F1_100x_check[which(F1_100x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F1.100x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F1.100x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F1.100x.opp),
                               row.names = "F1 Control - 100x")

# is there a more elegant way to do this?
for (i in colnames(all.df[2:16])){
  assign(i, data.frame(OTU = all.df$OTU,increasing = all.up[,i], decreasing = all.do[,i]))
}

## F1 Control - 10x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F1_10x = all.10x
opp_F1_10x_U <- opp_F1_10x[which(opp_F1_10x$F1_control_10x != 0 & opp_F1_10x$Unfiltered_control_10x != 0),]
opp_F1_10x_U <- opp_F1_10x_U[!(sign(opp_F1_10x_U$F1_control_10x) == sign(opp_F1_10x_U$Unfiltered_control_10x)),]

opp_F1_10x_U <- all.10x %>%
  subset(F1_control_10x != 0 & Unfiltered_control_10x != 0) %>%
  subset(!(sign(F1_control_10x) == sign(Unfiltered_control_10x)))

opp_F1_10x_F2 <- all.10x %>%
  subset(F1_control_10x != 0 & F2_control_10x != 0) %>%
  subset(!(sign(F1_control_10x) == sign(F2_control_10x))) 

opp_F1_10x_F3 <- all.10x %>%
  subset(F1_control_10x != 0 & F3_control_10x != 0) %>%
  subset(!(sign(F1_control_10x) == sign(F3_control_10x)))

opp_F1_10x_F4 <- all.10x %>%
  subset(F1_control_10x != 0 & F4_control_10x != 0) %>%
  subset(!(sign(F1_control_10x) == sign(F4_control_10x)))
opp_F1_10x_F4$OTU

all_opp_F1_10x = all.10x[,-2]
all_opp_F1_10x$Unfiltered_control_10x <- ifelse(all_opp_F1_10x$OTU %in% opp_F1_10x_U$OTU, 1, 0)
all_opp_F1_10x$F2_control_10x <- ifelse(all_opp_F1_10x$OTU %in% opp_F1_10x_F2$OTU, 1, 0)
all_opp_F1_10x$F3_control_10x <- ifelse(all_opp_F1_10x$OTU %in% opp_F1_10x_F3$OTU, 1, 0)
all_opp_F1_10x$F4_control_10x <- ifelse(all_opp_F1_10x$OTU %in% opp_F1_10x_F4$OTU, 1, 0)

diff_F1_10x = pres.nochange.10x
diff_F1_10x$Unfiltered_control_10x <- ifelse(diff_F1_10x$Unfiltered_control_10x == 1 | diff_F1_10x$OTU %in% opp_F1_10x_U$OTU, 1, 0)
diff_F1_10x$F2_control_10x <- ifelse(diff_F1_10x$F2_control_10x == 1 | diff_F1_10x$OTU %in% opp_F1_10x_F2$OTU, 1, 0)
diff_F1_10x$F3_control_10x <- ifelse(diff_F1_10x$F3_control_10x == 1 | diff_F1_10x$OTU %in% opp_F1_10x_F3$OTU, 1, 0)
diff_F1_10x$F4_control_10x <- ifelse(diff_F1_10x$F4_control_10x == 1 | diff_F1_10x$OTU %in% opp_F1_10x_F4$OTU, 1, 0)


F1_control_10x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F1_control_10x, decreasing = all.do$F1_control_10x,
                                    "Unfiltered Different" = diff_F1_10x$Unfiltered_control_10x,
                                    "F2 Different" = diff_F1_10x$F2_control_10x,
                                    "F3 Different" = diff_F1_10x$F3_control_10x,
                                    "F4 Different" = diff_F1_10x$F4_control_10x)
F1_10x_check <- subset(F1_control_10x_check, increasing + decreasing != 0)



nrow(F1_10x_check) #175
nrow(F1_10x_check[which(F1_10x_check$increasing == 1),]) #84
nrow(F1_10x_check[which(F1_10x_check$decreasing == 1),]) #91

# count how many direct
otus.F1.10x.up.direct <- F1_10x_check$OTU[which(rowSums(F1_10x_check[,-1]) == 1 & F1_10x_check$increasing == 1)] #7
otus.F1.10x.do.direct <- F1_10x_check$OTU[which(rowSums(F1_10x_check[,-1]) == 1 & F1_10x_check$decreasing == 1)] #10
otus.F1.10x.all.direct <- c(otus.F1.10x.up.direct, otus.F1.10x.do.direct) #17

#vector of indirect otus
`%!in%` <- Negate(`%in%`)
otus.F1.10x.indirect <- F1_10x_check$OTU[which(F1_10x_check$OTU %!in% otus.F1.10x.all.direct)] #158

# count opposite
otus.F1.10x.opp <- unique(c(opp_F1_10x_U$OTU, opp_F1_10x_F2$OTU, opp_F1_10x_F3$OTU, opp_F1_10x_F4$OTU)) #17

# put into data frame to add to affected.otus
affected.F1.10x <- data.frame("Indirect" = length(otus.F1.10x.indirect),
                              "OTUs Increasing" = nrow(F1_10x_check[which(F1_10x_check$increasing == 1),]),
                              "OTUs Decreasing" = nrow(F1_10x_check[which(F1_10x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F1.10x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F1.10x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F1.10x.opp),
                               row.names = "F1 Control - 10x")

affected.add <- rbind(affected.F1.100x, affected.F1.10x)

F1.10x.upset <- upset(F1_10x_check, sets = c("F4.Different", "F3.Different", "F2.Different", "Unfiltered.Different", "decreasing", "increasing"), 
        keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F1.10x.upset
grid.text(paste(length(F1_10x_check$OTU), "OTUs significantly affected in Fraction 1 10x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_10x_check[which(F1_10x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_10x_check[which(F1_10x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.10x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F1.10x$Increasing.in.all.fractions/affected.F1.10x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.10x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F1.10x$Decreasing.in.all.fractions/affected.F1.10x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F1.10x$Indirect, "of", length(F1_10x_check$OTU), "OTUs (", round(affected.F1.10x$Indirect/length(F1_10x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 


graph2ppt(x= F1.10x.upset ,file="out/F1.10x.upset.pptx")

## F1 Control - 1x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F1_1x = all.1x
opp_F1_1x_U <- opp_F1_1x[which(opp_F1_1x$F1_control_1x != 0 & opp_F1_1x$Unfiltered_control_1x != 0),]
opp_F1_1x_U <- opp_F1_1x_U[!(sign(opp_F1_1x_U$F1_control_1x) == sign(opp_F1_1x_U$Unfiltered_control_1x)),]

opp_F1_1x_U <- all.1x %>%
  subset(F1_control_1x != 0 & Unfiltered_control_1x != 0) %>%
  subset(!(sign(F1_control_1x) == sign(Unfiltered_control_1x)))

opp_F1_1x_F2 <- all.1x %>%
  subset(F1_control_1x != 0 & F2_control_1x != 0) %>%
  subset(!(sign(F1_control_1x) == sign(F2_control_1x))) 

opp_F1_1x_F3 <- all.1x %>%
  subset(F1_control_1x != 0 & F3_control_1x != 0) %>%
  subset(!(sign(F1_control_1x) == sign(F3_control_1x)))

opp_F1_1x_F4 <- all.1x %>%
  subset(F1_control_1x != 0 & F4_control_1x != 0) %>%
  subset(!(sign(F1_control_1x) == sign(F4_control_1x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F1_1x = pres.nochange.1x
diff_F1_1x$Unfiltered_control_1x <- ifelse(diff_F1_1x$Unfiltered_control_1x == 1 | diff_F1_1x$OTU %in% opp_F1_1x_U$OTU, 1, 0)
diff_F1_1x$F2_control_1x <- ifelse(diff_F1_1x$F2_control_1x == 1 | diff_F1_1x$OTU %in% opp_F1_1x_F2$OTU, 1, 0)
diff_F1_1x$F3_control_1x <- ifelse(diff_F1_1x$F3_control_1x == 1 | diff_F1_1x$OTU %in% opp_F1_1x_F3$OTU, 1, 0)
diff_F1_1x$F4_control_1x <- ifelse(diff_F1_1x$F4_control_1x == 1 | diff_F1_1x$OTU %in% opp_F1_1x_F4$OTU, 1, 0)

# create data frame for upset
F1_control_1x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F1_control_1x, decreasing = all.do$F1_control_1x,
                                   "Unfiltered Different" = diff_F1_1x$Unfiltered_control_1x,
                                   "F2 Different" = diff_F1_1x$F2_control_1x,
                                   "F3 Different" = diff_F1_1x$F3_control_1x,
                                   "F4 Different" = diff_F1_1x$F4_control_1x)
# get rid of OTUs that are not responding in this fraction
F1_1x_check <- subset(F1_control_1x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F1_1x_check) #60
nrow(F1_1x_check[which(F1_1x_check$increasing == 1),]) #21
nrow(F1_1x_check[which(F1_1x_check$decreasing == 1),]) #39

# vector of directly affected OTUs (the same response in all fractions)
otus.F1.1x.up.direct <- F1_1x_check$OTU[which(rowSums(F1_1x_check[,-1]) == 1 & F1_1x_check$increasing == 1)] #2
otus.F1.1x.do.direct <- F1_1x_check$OTU[which(rowSums(F1_1x_check[,-1]) == 1 & F1_1x_check$decreasing == 1)] #6
otus.F1.1x.all.direct <- c(otus.F1.1x.up.direct, otus.F1.1x.do.direct) #8

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F1.1x.indirect <- F1_1x_check$OTU[which(F1_1x_check$OTU %!in% otus.F1.1x.all.direct)] #52

# count opposite
otus.F1.1x.opp <- unique(c(opp_F1_1x_U$OTU, opp_F1_1x_F2$OTU, opp_F1_1x_F3$OTU, opp_F1_1x_F4$OTU)) #6

# put into data frame to add to affected.otus
affected.F1.1x <- data.frame("Indirect" = length(otus.F1.1x.indirect),
                              "OTUs Increasing" = nrow(F1_1x_check[which(F1_1x_check$increasing == 1),]),
                              "OTUs Decreasing" = nrow(F1_1x_check[which(F1_1x_check$decreasing == 1),]),
                              "Increasing in all fractions" = length(otus.F1.1x.up.direct), 
                              "Decreasing in all fractions" = length(otus.F1.1x.do.direct),
                              "Opposite in at least one fraction" = length(otus.F1.1x.opp),
                              row.names = "F1 Control - 1x")
affected.add <- rbind(affected.add, affected.F1.1x)

# Make upset plot
F1.1x.upset <- upset(F1_1x_check, sets = c("F4.Different", "F3.Different", "F2.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                      keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F1.1x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F1_1x_check$OTU), "OTUs significantly affected in Fraction 1 1x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_1x_check[which(F1_1x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F1_1x_check[which(F1_1x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.1x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F1.1x$Increasing.in.all.fractions/affected.F1.1x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F1.1x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F1.1x$Decreasing.in.all.fractions/affected.F1.1x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F1.1x$Indirect, "of", length(F1_1x_check$OTU), "OTUs (", round(affected.F1.1x$Indirect/length(F1_1x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F1.1x.upset ,file="out/F1.1x.upset.pptx")

## F2 Control - 100x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F2_100x = all.100x
opp_F2_100x_U <- opp_F2_100x[which(opp_F2_100x$F2_control_100x != 0 & opp_F2_100x$Unfiltered_control_100x != 0),]
opp_F2_100x_U <- opp_F2_100x_U[!(sign(opp_F2_100x_U$F2_control_100x) == sign(opp_F2_100x_U$Unfiltered_control_100x)),]

opp_F2_100x_U <- all.100x %>%
  subset(F2_control_100x != 0 & Unfiltered_control_100x != 0) %>%
  subset(!(sign(F2_control_100x) == sign(Unfiltered_control_100x)))

opp_F2_100x_F1 <- all.100x %>%
  subset(F2_control_100x != 0 & F1_control_100x != 0) %>%
  subset(!(sign(F2_control_100x) == sign(F1_control_100x))) 

opp_F2_100x_F3 <- all.100x %>%
  subset(F2_control_100x != 0 & F3_control_100x != 0) %>%
  subset(!(sign(F2_control_100x) == sign(F3_control_100x)))

opp_F2_100x_F4 <- all.100x %>%
  subset(F2_control_100x != 0 & F4_control_100x != 0) %>%
  subset(!(sign(F2_control_100x) == sign(F4_control_100x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F2_100x = pres.nochange.100x
diff_F2_100x$Unfiltered_control_100x <- ifelse(diff_F2_100x$Unfiltered_control_100x == 1 | diff_F2_100x$OTU %in% opp_F2_100x_U$OTU, 1, 0)
diff_F2_100x$F1_control_100x <- ifelse(diff_F2_100x$F1_control_100x == 1 | diff_F2_100x$OTU %in% opp_F2_100x_F1$OTU, 1, 0)
diff_F2_100x$F3_control_100x <- ifelse(diff_F2_100x$F3_control_100x == 1 | diff_F2_100x$OTU %in% opp_F2_100x_F3$OTU, 1, 0)
diff_F2_100x$F4_control_100x <- ifelse(diff_F2_100x$F4_control_100x == 1 | diff_F2_100x$OTU %in% opp_F2_100x_F4$OTU, 1, 0)

# create data frame for upset
F2_control_100x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F2_control_100x, decreasing = all.do$F2_control_100x,
                                  "Unfiltered Different" = diff_F2_100x$Unfiltered_control_100x,
                                  "F1 Different" = diff_F2_100x$F1_control_100x,
                                  "F3 Different" = diff_F2_100x$F3_control_100x,
                                  "F4 Different" = diff_F2_100x$F4_control_100x)
# get rid of OTUs that are not responding in this fraction
F2_100x_check <- subset(F2_control_100x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F2_100x_check) #133
nrow(F2_100x_check[which(F2_100x_check$increasing == 1),]) #65
nrow(F2_100x_check[which(F2_100x_check$decreasing == 1),]) #68

# vector of directly affected OTUs (the same response in all fractions)
otus.F2.100x.up.direct <- F2_100x_check$OTU[which(rowSums(F2_100x_check[,-1]) == 1 & F2_100x_check$increasing == 1)] #35
otus.F2.100x.do.direct <- F2_100x_check$OTU[which(rowSums(F2_100x_check[,-1]) == 1 & F2_100x_check$decreasing == 1)] #20
otus.F2.100x.all.direct <- c(otus.F2.100x.up.direct, otus.F2.100x.do.direct) #55

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F2.100x.indirect <- F2_100x_check$OTU[which(F2_100x_check$OTU %!in% otus.F2.100x.all.direct)] #78

# count opposite
otus.F2.100x.opp <- unique(c(opp_F2_100x_U$OTU, opp_F2_100x_F1$OTU, opp_F2_100x_F3$OTU, opp_F2_100x_F4$OTU)) #20

# put into data frame to add to affected.otus
affected.F2.100x <- data.frame("Indirect" = length(otus.F2.100x.indirect),
                             "OTUs Increasing" = nrow(F2_100x_check[which(F2_100x_check$increasing == 1),]),
                             "OTUs Decreasing" = nrow(F2_100x_check[which(F2_100x_check$decreasing == 1),]),
                             "Increasing in all fractions" = length(otus.F2.100x.up.direct), 
                             "Decreasing in all fractions" = length(otus.F2.100x.do.direct),
                             "Opposite in at least one fraction" = length(otus.F2.100x.opp),
                             row.names = "F2 Control - 100x")
affected.add <- rbind(affected.add, affected.F2.100x)

# Make upset plot
F2.100x.upset <- upset(F2_100x_check, sets = c("F4.Different", "F3.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                     keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F2.100x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F2_100x_check$OTU), "OTUs significantly affected in F2 100x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_100x_check[which(F2_100x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_100x_check[which(F2_100x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.100x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F2.100x$Increasing.in.all.fractions/affected.F2.100x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.100x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F2.100x$Decreasing.in.all.fractions/affected.F2.100x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F2.100x$Indirect, "of", length(F2_100x_check$OTU), "OTUs (", round(affected.F2.100x$Indirect/length(F2_100x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F2.100x.upset ,file="out/F2.100x.upset.pptx")

## F2 Control - 10x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F2_10x = all.10x
opp_F2_10x_U <- opp_F2_10x[which(opp_F2_10x$F2_control_10x != 0 & opp_F2_10x$Unfiltered_control_10x != 0),]
opp_F2_10x_U <- opp_F2_10x_U[!(sign(opp_F2_10x_U$F2_control_10x) == sign(opp_F2_10x_U$Unfiltered_control_10x)),]

opp_F2_10x_U <- all.10x %>%
  subset(F2_control_10x != 0 & Unfiltered_control_10x != 0) %>%
  subset(!(sign(F2_control_10x) == sign(Unfiltered_control_10x)))

opp_F2_10x_F1 <- all.10x %>%
  subset(F2_control_10x != 0 & F1_control_10x != 0) %>%
  subset(!(sign(F2_control_10x) == sign(F1_control_10x))) 

opp_F2_10x_F3 <- all.10x %>%
  subset(F2_control_10x != 0 & F3_control_10x != 0) %>%
  subset(!(sign(F2_control_10x) == sign(F3_control_10x)))

opp_F2_10x_F4 <- all.10x %>%
  subset(F2_control_10x != 0 & F4_control_10x != 0) %>%
  subset(!(sign(F2_control_10x) == sign(F4_control_10x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F2_10x = pres.nochange.10x
diff_F2_10x$Unfiltered_control_10x <- ifelse(diff_F2_10x$Unfiltered_control_10x == 1 | diff_F2_10x$OTU %in% opp_F2_10x_U$OTU, 1, 0)
diff_F2_10x$F1_control_10x <- ifelse(diff_F2_10x$F1_control_10x == 1 | diff_F2_10x$OTU %in% opp_F2_10x_F1$OTU, 1, 0)
diff_F2_10x$F3_control_10x <- ifelse(diff_F2_10x$F3_control_10x == 1 | diff_F2_10x$OTU %in% opp_F2_10x_F3$OTU, 1, 0)
diff_F2_10x$F4_control_10x <- ifelse(diff_F2_10x$F4_control_10x == 1 | diff_F2_10x$OTU %in% opp_F2_10x_F4$OTU, 1, 0)

# create data frame for upset
F2_control_10x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F2_control_10x, decreasing = all.do$F2_control_10x,
                                    "Unfiltered Different" = diff_F2_10x$Unfiltered_control_10x,
                                    "F1 Different" = diff_F2_10x$F1_control_10x,
                                    "F3 Different" = diff_F2_10x$F3_control_10x,
                                    "F4 Different" = diff_F2_10x$F4_control_10x)
# get rid of OTUs that are not responding in this fraction
F2_10x_check <- subset(F2_control_10x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F2_10x_check) #91
nrow(F2_10x_check[which(F2_10x_check$increasing == 1),]) #39
nrow(F2_10x_check[which(F2_10x_check$decreasing == 1),]) #52

# vector of directly affected OTUs (the same response in all fractions)
otus.F2.10x.up.direct <- F2_10x_check$OTU[which(rowSums(F2_10x_check[,-1]) == 1 & F2_10x_check$increasing == 1)] #7
otus.F2.10x.do.direct <- F2_10x_check$OTU[which(rowSums(F2_10x_check[,-1]) == 1 & F2_10x_check$decreasing == 1)] #9
otus.F2.10x.all.direct <- c(otus.F2.10x.up.direct, otus.F2.10x.do.direct) #16

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F2.10x.indirect <- F2_10x_check$OTU[which(F2_10x_check$OTU %!in% otus.F2.10x.all.direct)] #75

# count opposite
otus.F2.10x.opp <- unique(c(opp_F2_10x_U$OTU, opp_F2_10x_F1$OTU, opp_F2_10x_F3$OTU, opp_F2_10x_F4$OTU)) #16

# put into data frame to add to affected.otus
affected.F2.10x <- data.frame("Indirect" = length(otus.F2.10x.indirect),
                               "OTUs Increasing" = nrow(F2_10x_check[which(F2_10x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(F2_10x_check[which(F2_10x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F2.10x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F2.10x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F2.10x.opp),
                               row.names = "F2 Control - 10x")
affected.add <- rbind(affected.add, affected.F2.10x)

# Make upset plot
F2.10x.upset <- upset(F2_10x_check, sets = c("F4.Different", "F3.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                       keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F2.10x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F2_10x_check$OTU), "OTUs significantly affected in F2 10x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_10x_check[which(F2_10x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_10x_check[which(F2_10x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.10x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F2.10x$Increasing.in.all.fractions/affected.F2.10x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.10x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F2.10x$Decreasing.in.all.fractions/affected.F2.10x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F2.10x$Indirect, "of", length(F2_10x_check$OTU), "OTUs (", round(affected.F2.10x$Indirect/length(F2_10x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F2.10x.upset ,file="out/F2.10x.upset.pptx")

## F2 Control - 1x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F2_1x = all.1x
opp_F2_1x_U <- opp_F2_1x[which(opp_F2_1x$F2_control_1x != 0 & opp_F2_1x$Unfiltered_control_1x != 0),]
opp_F2_1x_U <- opp_F2_1x_U[!(sign(opp_F2_1x_U$F2_control_1x) == sign(opp_F2_1x_U$Unfiltered_control_1x)),]

opp_F2_1x_U <- all.1x %>%
  subset(F2_control_1x != 0 & Unfiltered_control_1x != 0) %>%
  subset(!(sign(F2_control_1x) == sign(Unfiltered_control_1x)))

opp_F2_1x_F1 <- all.1x %>%
  subset(F2_control_1x != 0 & F1_control_1x != 0) %>%
  subset(!(sign(F2_control_1x) == sign(F1_control_1x))) 

opp_F2_1x_F3 <- all.1x %>%
  subset(F2_control_1x != 0 & F3_control_1x != 0) %>%
  subset(!(sign(F2_control_1x) == sign(F3_control_1x)))

opp_F2_1x_F4 <- all.1x %>%
  subset(F2_control_1x != 0 & F4_control_1x != 0) %>%
  subset(!(sign(F2_control_1x) == sign(F4_control_1x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F2_1x = pres.nochange.1x
diff_F2_1x$Unfiltered_control_1x <- ifelse(diff_F2_1x$Unfiltered_control_1x == 1 | diff_F2_1x$OTU %in% opp_F2_1x_U$OTU, 1, 0)
diff_F2_1x$F1_control_1x <- ifelse(diff_F2_1x$F1_control_1x == 1 | diff_F2_1x$OTU %in% opp_F2_1x_F1$OTU, 1, 0)
diff_F2_1x$F3_control_1x <- ifelse(diff_F2_1x$F3_control_1x == 1 | diff_F2_1x$OTU %in% opp_F2_1x_F3$OTU, 1, 0)
diff_F2_1x$F4_control_1x <- ifelse(diff_F2_1x$F4_control_1x == 1 | diff_F2_1x$OTU %in% opp_F2_1x_F4$OTU, 1, 0)

# create data frame for upset
F2_control_1x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F2_control_1x, decreasing = all.do$F2_control_1x,
                                    "Unfiltered Different" = diff_F2_1x$Unfiltered_control_1x,
                                    "F1 Different" = diff_F2_1x$F1_control_1x,
                                    "F3 Different" = diff_F2_1x$F3_control_1x,
                                    "F4 Different" = diff_F2_1x$F4_control_1x)
# get rid of OTUs that are not responding in this fraction
F2_1x_check <- subset(F2_control_1x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F2_1x_check) #56
nrow(F2_1x_check[which(F2_1x_check$increasing == 1),]) #25
nrow(F2_1x_check[which(F2_1x_check$decreasing == 1),]) #31

# vector of directly affected OTUs (the same response in all fractions)
otus.F2.1x.up.direct <- F2_1x_check$OTU[which(rowSums(F2_1x_check[,-1]) == 1 & F2_1x_check$increasing == 1)] #2
otus.F2.1x.do.direct <- F2_1x_check$OTU[which(rowSums(F2_1x_check[,-1]) == 1 & F2_1x_check$decreasing == 1)] #5
otus.F2.1x.all.direct <- c(otus.F2.1x.up.direct, otus.F2.1x.do.direct) #7

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F2.1x.indirect <- F2_1x_check$OTU[which(F2_1x_check$OTU %!in% otus.F2.1x.all.direct)] #49

# count opposite
otus.F2.1x.opp <- unique(c(opp_F2_1x_U$OTU, opp_F2_1x_F1$OTU, opp_F2_1x_F3$OTU, opp_F2_1x_F4$OTU)) #5

# put into data frame to add to affected.otus
affected.F2.1x <- data.frame("Indirect" = length(otus.F2.1x.indirect),
                               "OTUs Increasing" = nrow(F2_1x_check[which(F2_1x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(F2_1x_check[which(F2_1x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F2.1x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F2.1x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F2.1x.opp),
                               row.names = "F2 Control - 1x")
affected.add <- rbind(affected.add, affected.F2.1x)

# Make upset plot
F2.1x.upset <- upset(F2_1x_check, sets = c("F4.Different", "F3.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                       keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F2.1x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F2_1x_check$OTU), "OTUs significantly affected in F2 1x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_1x_check[which(F2_1x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F2_1x_check[which(F2_1x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.1x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F2.1x$Increasing.in.all.fractions/affected.F2.1x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F2.1x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F2.1x$Decreasing.in.all.fractions/affected.F2.1x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F2.1x$Indirect, "of", length(F2_1x_check$OTU), "OTUs (", round(affected.F2.1x$Indirect/length(F2_1x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F2.1x.upset ,file="out/F2.1x.upset.pptx")

## F3 Control - 100x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F3_100x = all.100x
opp_F3_100x_U <- opp_F3_100x[which(opp_F3_100x$F3_control_100x != 0 & opp_F3_100x$Unfiltered_control_100x != 0),]
opp_F3_100x_U <- opp_F3_100x_U[!(sign(opp_F3_100x_U$F3_control_100x) == sign(opp_F3_100x_U$Unfiltered_control_100x)),]

opp_F3_100x_U <- all.100x %>%
  subset(F3_control_100x != 0 & Unfiltered_control_100x != 0) %>%
  subset(!(sign(F3_control_100x) == sign(Unfiltered_control_100x)))

opp_F3_100x_F1 <- all.100x %>%
  subset(F3_control_100x != 0 & F1_control_100x != 0) %>%
  subset(!(sign(F3_control_100x) == sign(F1_control_100x))) 

opp_F3_100x_F2 <- all.100x %>%
  subset(F3_control_100x != 0 & F2_control_100x != 0) %>%
  subset(!(sign(F3_control_100x) == sign(F2_control_100x)))

opp_F3_100x_F4 <- all.100x %>%
  subset(F3_control_100x != 0 & F4_control_100x != 0) %>%
  subset(!(sign(F3_control_100x) == sign(F4_control_100x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F3_100x = pres.nochange.100x
diff_F3_100x$Unfiltered_control_100x <- ifelse(diff_F3_100x$Unfiltered_control_100x == 1 | diff_F3_100x$OTU %in% opp_F3_100x_U$OTU, 1, 0)
diff_F3_100x$F1_control_100x <- ifelse(diff_F3_100x$F1_control_100x == 1 | diff_F3_100x$OTU %in% opp_F3_100x_F1$OTU, 1, 0)
diff_F3_100x$F2_control_100x <- ifelse(diff_F3_100x$F2_control_100x == 1 | diff_F3_100x$OTU %in% opp_F3_100x_F2$OTU, 1, 0)
diff_F3_100x$F4_control_100x <- ifelse(diff_F3_100x$F4_control_100x == 1 | diff_F3_100x$OTU %in% opp_F3_100x_F4$OTU, 1, 0)

# create data frame for upset
F3_control_100x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F3_control_100x, decreasing = all.do$F3_control_100x,
                                    "Unfiltered Different" = diff_F3_100x$Unfiltered_control_100x,
                                    "F1 Different" = diff_F3_100x$F1_control_100x,
                                    "F2 Different" = diff_F3_100x$F2_control_100x,
                                    "F4 Different" = diff_F3_100x$F4_control_100x)
# get rid of OTUs that are not responding in this fraction
F3_100x_check <- subset(F3_control_100x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F3_100x_check) #209
nrow(F3_100x_check[which(F3_100x_check$increasing == 1),]) #104
nrow(F3_100x_check[which(F3_100x_check$decreasing == 1),]) #105

# vector of directly affected OTUs (the same response in all fractions)
otus.F3.100x.up.direct <- F3_100x_check$OTU[which(rowSums(F3_100x_check[,-1]) == 1 & F3_100x_check$increasing == 1)] #40
otus.F3.100x.do.direct <- F3_100x_check$OTU[which(rowSums(F3_100x_check[,-1]) == 1 & F3_100x_check$decreasing == 1)] #22
otus.F3.100x.all.direct <- c(otus.F3.100x.up.direct, otus.F3.100x.do.direct) #62

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F3.100x.indirect <- F3_100x_check$OTU[which(F3_100x_check$OTU %!in% otus.F3.100x.all.direct)] #147

# count opposite
otus.F3.100x.opp <- unique(c(opp_F3_100x_U$OTU, opp_F3_100x_F1$OTU, opp_F3_100x_F2$OTU, opp_F3_100x_F4$OTU)) #22

# put into data frame to add to affected.otus
affected.F3.100x <- data.frame("Indirect" = length(otus.F3.100x.indirect),
                               "OTUs Increasing" = nrow(F3_100x_check[which(F3_100x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(F3_100x_check[which(F3_100x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F3.100x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F3.100x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F3.100x.opp),
                               row.names = "F3 Control - 100x")
affected.add <- rbind(affected.add, affected.F3.100x)

# Make upset plot
F3.100x.upset <- upset(F3_100x_check, sets = c("F4.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                       keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F3.100x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F3_100x_check$OTU), "OTUs significantly affected in F3 100x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_100x_check[which(F3_100x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_100x_check[which(F3_100x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.100x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F3.100x$Increasing.in.all.fractions/affected.F3.100x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.100x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F3.100x$Decreasing.in.all.fractions/affected.F3.100x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F3.100x$Indirect, "of", length(F3_100x_check$OTU), "OTUs (", round(affected.F3.100x$Indirect/length(F3_100x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F3.100x.upset ,file="out/F3.100x.upset.pptx")

## F3 Control - 10x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F3_10x = all.10x
opp_F3_10x_U <- opp_F3_10x[which(opp_F3_10x$F3_control_10x != 0 & opp_F3_10x$Unfiltered_control_10x != 0),]
opp_F3_10x_U <- opp_F3_10x_U[!(sign(opp_F3_10x_U$F3_control_10x) == sign(opp_F3_10x_U$Unfiltered_control_10x)),]

opp_F3_10x_U <- all.10x %>%
  subset(F3_control_10x != 0 & Unfiltered_control_10x != 0) %>%
  subset(!(sign(F3_control_10x) == sign(Unfiltered_control_10x)))

opp_F3_10x_F1 <- all.10x %>%
  subset(F3_control_10x != 0 & F1_control_10x != 0) %>%
  subset(!(sign(F3_control_10x) == sign(F1_control_10x))) 

opp_F3_10x_F2 <- all.10x %>%
  subset(F3_control_10x != 0 & F2_control_10x != 0) %>%
  subset(!(sign(F3_control_10x) == sign(F2_control_10x)))

opp_F3_10x_F4 <- all.10x %>%
  subset(F3_control_10x != 0 & F4_control_10x != 0) %>%
  subset(!(sign(F3_control_10x) == sign(F4_control_10x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F3_10x = pres.nochange.10x
diff_F3_10x$Unfiltered_control_10x <- ifelse(diff_F3_10x$Unfiltered_control_10x == 1 | diff_F3_10x$OTU %in% opp_F3_10x_U$OTU, 1, 0)
diff_F3_10x$F1_control_10x <- ifelse(diff_F3_10x$F1_control_10x == 1 | diff_F3_10x$OTU %in% opp_F3_10x_F1$OTU, 1, 0)
diff_F3_10x$F2_control_10x <- ifelse(diff_F3_10x$F2_control_10x == 1 | diff_F3_10x$OTU %in% opp_F3_10x_F2$OTU, 1, 0)
diff_F3_10x$F4_control_10x <- ifelse(diff_F3_10x$F4_control_10x == 1 | diff_F3_10x$OTU %in% opp_F3_10x_F4$OTU, 1, 0)

# create data frame for upset
F3_control_10x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F3_control_10x, decreasing = all.do$F3_control_10x,
                                   "Unfiltered Different" = diff_F3_10x$Unfiltered_control_10x,
                                   "F1 Different" = diff_F3_10x$F1_control_10x,
                                   "F2 Different" = diff_F3_10x$F2_control_10x,
                                   "F4 Different" = diff_F3_10x$F4_control_10x)
# get rid of OTUs that are not responding in this fraction
F3_10x_check <- subset(F3_control_10x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F3_10x_check) #165
nrow(F3_10x_check[which(F3_10x_check$increasing == 1),]) #65
nrow(F3_10x_check[which(F3_10x_check$decreasing == 1),]) #100

# vector of directly affected OTUs (the same response in all fractions)
otus.F3.10x.up.direct <- F3_10x_check$OTU[which(rowSums(F3_10x_check[,-1]) == 1 & F3_10x_check$increasing == 1)] #7
otus.F3.10x.do.direct <- F3_10x_check$OTU[which(rowSums(F3_10x_check[,-1]) == 1 & F3_10x_check$decreasing == 1)] #10
otus.F3.10x.all.direct <- c(otus.F3.10x.up.direct, otus.F3.10x.do.direct) #17

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F3.10x.indirect <- F3_10x_check$OTU[which(F3_10x_check$OTU %!in% otus.F3.10x.all.direct)] #148

# count opposite
otus.F3.10x.opp <- unique(c(opp_F3_10x_U$OTU, opp_F3_10x_F1$OTU, opp_F3_10x_F2$OTU, opp_F3_10x_F4$OTU)) #22

# put into data frame to add to affected.otus
affected.F3.10x <- data.frame("Indirect" = length(otus.F3.10x.indirect),
                              "OTUs Increasing" = nrow(F3_10x_check[which(F3_10x_check$increasing == 1),]),
                              "OTUs Decreasing" = nrow(F3_10x_check[which(F3_10x_check$decreasing == 1),]),
                              "Increasing in all fractions" = length(otus.F3.10x.up.direct), 
                              "Decreasing in all fractions" = length(otus.F3.10x.do.direct),
                              "Opposite in at least one fraction" = length(otus.F3.10x.opp),
                              row.names = "F3 Control - 10x")
affected.add <- rbind(affected.add, affected.F3.10x)

# Make upset plot
F3.10x.upset <- upset(F3_10x_check, sets = c("F4.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                      keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F3.10x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F3_10x_check$OTU), "OTUs significantly affected in F3 10x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_10x_check[which(F3_10x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_10x_check[which(F3_10x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.10x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F3.10x$Increasing.in.all.fractions/affected.F3.10x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.10x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F3.10x$Decreasing.in.all.fractions/affected.F3.10x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F3.10x$Indirect, "of", length(F3_10x_check$OTU), "OTUs (", round(affected.F3.10x$Indirect/length(F3_10x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F3.10x.upset ,file="out/F3.10x.upset.pptx")

## F3 Control - 1x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F3_1x = all.1x
opp_F3_1x_U <- opp_F3_1x[which(opp_F3_1x$F3_control_1x != 0 & opp_F3_1x$Unfiltered_control_1x != 0),]
opp_F3_1x_U <- opp_F3_1x_U[!(sign(opp_F3_1x_U$F3_control_1x) == sign(opp_F3_1x_U$Unfiltered_control_1x)),]

opp_F3_1x_U <- all.1x %>%
  subset(F3_control_1x != 0 & Unfiltered_control_1x != 0) %>%
  subset(!(sign(F3_control_1x) == sign(Unfiltered_control_1x)))

opp_F3_1x_F1 <- all.1x %>%
  subset(F3_control_1x != 0 & F1_control_1x != 0) %>%
  subset(!(sign(F3_control_1x) == sign(F1_control_1x))) 

opp_F3_1x_F2 <- all.1x %>%
  subset(F3_control_1x != 0 & F2_control_1x != 0) %>%
  subset(!(sign(F3_control_1x) == sign(F2_control_1x)))

opp_F3_1x_F4 <- all.1x %>%
  subset(F3_control_1x != 0 & F4_control_1x != 0) %>%
  subset(!(sign(F3_control_1x) == sign(F4_control_1x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F3_1x = pres.nochange.1x
diff_F3_1x$Unfiltered_control_1x <- ifelse(diff_F3_1x$Unfiltered_control_1x == 1 | diff_F3_1x$OTU %in% opp_F3_1x_U$OTU, 1, 0)
diff_F3_1x$F1_control_1x <- ifelse(diff_F3_1x$F1_control_1x == 1 | diff_F3_1x$OTU %in% opp_F3_1x_F1$OTU, 1, 0)
diff_F3_1x$F2_control_1x <- ifelse(diff_F3_1x$F2_control_1x == 1 | diff_F3_1x$OTU %in% opp_F3_1x_F2$OTU, 1, 0)
diff_F3_1x$F4_control_1x <- ifelse(diff_F3_1x$F4_control_1x == 1 | diff_F3_1x$OTU %in% opp_F3_1x_F4$OTU, 1, 0)

# create data frame for upset
F3_control_1x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F3_control_1x, decreasing = all.do$F3_control_1x,
                                  "Unfiltered Different" = diff_F3_1x$Unfiltered_control_1x,
                                  "F1 Different" = diff_F3_1x$F1_control_1x,
                                  "F2 Different" = diff_F3_1x$F2_control_1x,
                                  "F4 Different" = diff_F3_1x$F4_control_1x)
# get rid of OTUs that are not responding in this fraction
F3_1x_check <- subset(F3_control_1x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F3_1x_check) #102
nrow(F3_1x_check[which(F3_1x_check$increasing == 1),]) #32
nrow(F3_1x_check[which(F3_1x_check$decreasing == 1),]) #70

# vector of directly affected OTUs (the same response in all fractions)
otus.F3.1x.up.direct <- F3_1x_check$OTU[which(rowSums(F3_1x_check[,-1]) == 1 & F3_1x_check$increasing == 1)] #2
otus.F3.1x.do.direct <- F3_1x_check$OTU[which(rowSums(F3_1x_check[,-1]) == 1 & F3_1x_check$decreasing == 1)] #6
otus.F3.1x.all.direct <- c(otus.F3.1x.up.direct, otus.F3.1x.do.direct) #8

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F3.1x.indirect <- F3_1x_check$OTU[which(F3_1x_check$OTU %!in% otus.F3.1x.all.direct)] #94

# count opposite
otus.F3.1x.opp <- unique(c(opp_F3_1x_U$OTU, opp_F3_1x_F1$OTU, opp_F3_1x_F2$OTU, opp_F3_1x_F4$OTU)) #15

# put into data frame to add to affected.otus
affected.F3.1x <- data.frame("Indirect" = length(otus.F3.1x.indirect),
                             "OTUs Increasing" = nrow(F3_1x_check[which(F3_1x_check$increasing == 1),]),
                             "OTUs Decreasing" = nrow(F3_1x_check[which(F3_1x_check$decreasing == 1),]),
                             "Increasing in all fractions" = length(otus.F3.1x.up.direct), 
                             "Decreasing in all fractions" = length(otus.F3.1x.do.direct),
                             "Opposite in at least one fraction" = length(otus.F3.1x.opp),
                             row.names = "F3 Control - 1x")
affected.add <- rbind(affected.add, affected.F3.1x)

# Make upset plot
F3.1x.upset <- upset(F3_1x_check, sets = c("F4.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                     keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F3.1x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F3_1x_check$OTU), "OTUs significantly affected in F3 1x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_1x_check[which(F3_1x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F3_1x_check[which(F3_1x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.1x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F3.1x$Increasing.in.all.fractions/affected.F3.1x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F3.1x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F3.1x$Decreasing.in.all.fractions/affected.F3.1x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F3.1x$Indirect, "of", length(F3_1x_check$OTU), "OTUs (", round(affected.F3.1x$Indirect/length(F3_1x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F3.1x.upset ,file="out/F3.1x.upset.pptx")

## F4 Control - 100x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F4_100x = all.100x
opp_F4_100x_U <- opp_F4_100x[which(opp_F4_100x$F4_control_100x != 0 & opp_F4_100x$Unfiltered_control_100x != 0),]
opp_F4_100x_U <- opp_F4_100x_U[!(sign(opp_F4_100x_U$F4_control_100x) == sign(opp_F4_100x_U$Unfiltered_control_100x)),]

opp_F4_100x_U <- all.100x %>%
  subset(F4_control_100x != 0 & Unfiltered_control_100x != 0) %>%
  subset(!(sign(F4_control_100x) == sign(Unfiltered_control_100x)))

opp_F4_100x_F1 <- all.100x %>%
  subset(F4_control_100x != 0 & F1_control_100x != 0) %>%
  subset(!(sign(F4_control_100x) == sign(F1_control_100x))) 

opp_F4_100x_F2 <- all.100x %>%
  subset(F4_control_100x != 0 & F2_control_100x != 0) %>%
  subset(!(sign(F4_control_100x) == sign(F2_control_100x)))

opp_F4_100x_F3 <- all.100x %>%
  subset(F4_control_100x != 0 & F3_control_100x != 0) %>%
  subset(!(sign(F4_control_100x) == sign(F3_control_100x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F4_100x = pres.nochange.100x
diff_F4_100x$Unfiltered_control_100x <- ifelse(diff_F4_100x$Unfiltered_control_100x == 1 | diff_F4_100x$OTU %in% opp_F4_100x_U$OTU, 1, 0)
diff_F4_100x$F1_control_100x <- ifelse(diff_F4_100x$F1_control_100x == 1 | diff_F4_100x$OTU %in% opp_F4_100x_F1$OTU, 1, 0)
diff_F4_100x$F2_control_100x <- ifelse(diff_F4_100x$F2_control_100x == 1 | diff_F4_100x$OTU %in% opp_F4_100x_F2$OTU, 1, 0)
diff_F4_100x$F3_control_100x <- ifelse(diff_F4_100x$F3_control_100x == 1 | diff_F4_100x$OTU %in% opp_F4_100x_F3$OTU, 1, 0)

# create data frame for upset
F4_control_100x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F4_control_100x, decreasing = all.do$F4_control_100x,
                                    "Unfiltered Different" = diff_F4_100x$Unfiltered_control_100x,
                                    "F1 Different" = diff_F4_100x$F1_control_100x,
                                    "F2 Different" = diff_F4_100x$F2_control_100x,
                                    "F3 Different" = diff_F4_100x$F3_control_100x)
# get rid of OTUs that are not responding in this fraction
F4_100x_check <- subset(F4_control_100x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F4_100x_check) #143
nrow(F4_100x_check[which(F4_100x_check$increasing == 1),]) #65
nrow(F4_100x_check[which(F4_100x_check$decreasing == 1),]) #78

# vector of directly affected OTUs (the same response in all fractions)
otus.F4.100x.up.direct <- F4_100x_check$OTU[which(rowSums(F4_100x_check[,-1]) == 1 & F4_100x_check$increasing == 1)] #15
otus.F4.100x.do.direct <- F4_100x_check$OTU[which(rowSums(F4_100x_check[,-1]) == 1 & F4_100x_check$decreasing == 1)] #20
otus.F4.100x.all.direct <- c(otus.F4.100x.up.direct, otus.F4.100x.do.direct) #35

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F4.100x.indirect <- F4_100x_check$OTU[which(F4_100x_check$OTU %!in% otus.F4.100x.all.direct)] #108

# count opposite
otus.F4.100x.opp <- unique(c(opp_F4_100x_U$OTU, opp_F4_100x_F1$OTU, opp_F4_100x_F2$OTU, opp_F4_100x_F3$OTU)) #21

# put into data frame to add to affected.otus
affected.F4.100x <- data.frame("Indirect" = length(otus.F4.100x.indirect),
                               "OTUs Increasing" = nrow(F4_100x_check[which(F4_100x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(F4_100x_check[which(F4_100x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.F4.100x.up.direct), 
                               "Decreasing in all fractions" = length(otus.F4.100x.do.direct),
                               "Opposite in at least one fraction" = length(otus.F4.100x.opp),
                               row.names = "F4 Control - 100x")
affected.add <- rbind(affected.add, affected.F4.100x)

# Make upset plot
F4.100x.upset <- upset(F4_100x_check, sets = c("F3.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                       keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F4.100x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F4_100x_check$OTU), "OTUs significantly affected in F4 100x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_100x_check[which(F4_100x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_100x_check[which(F4_100x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.100x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F4.100x$Increasing.in.all.fractions/affected.F4.100x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.100x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F4.100x$Decreasing.in.all.fractions/affected.F4.100x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F4.100x$Indirect, "of", length(F4_100x_check$OTU), "OTUs (", round(affected.F4.100x$Indirect/length(F4_100x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F4.100x.upset ,file="out/F4.100x.upset.pptx")

## F4 Control - 10x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F4_10x = all.10x
opp_F4_10x_U <- opp_F4_10x[which(opp_F4_10x$F4_control_10x != 0 & opp_F4_10x$Unfiltered_control_10x != 0),]
opp_F4_10x_U <- opp_F4_10x_U[!(sign(opp_F4_10x_U$F4_control_10x) == sign(opp_F4_10x_U$Unfiltered_control_10x)),]

opp_F4_10x_U <- all.10x %>%
  subset(F4_control_10x != 0 & Unfiltered_control_10x != 0) %>%
  subset(!(sign(F4_control_10x) == sign(Unfiltered_control_10x)))

opp_F4_10x_F1 <- all.10x %>%
  subset(F4_control_10x != 0 & F1_control_10x != 0) %>%
  subset(!(sign(F4_control_10x) == sign(F1_control_10x))) 

opp_F4_10x_F2 <- all.10x %>%
  subset(F4_control_10x != 0 & F2_control_10x != 0) %>%
  subset(!(sign(F4_control_10x) == sign(F2_control_10x)))

opp_F4_10x_F3 <- all.10x %>%
  subset(F4_control_10x != 0 & F3_control_10x != 0) %>%
  subset(!(sign(F4_control_10x) == sign(F3_control_10x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F4_10x = pres.nochange.10x
diff_F4_10x$Unfiltered_control_10x <- ifelse(diff_F4_10x$Unfiltered_control_10x == 1 | diff_F4_10x$OTU %in% opp_F4_10x_U$OTU, 1, 0)
diff_F4_10x$F1_control_10x <- ifelse(diff_F4_10x$F1_control_10x == 1 | diff_F4_10x$OTU %in% opp_F4_10x_F1$OTU, 1, 0)
diff_F4_10x$F2_control_10x <- ifelse(diff_F4_10x$F2_control_10x == 1 | diff_F4_10x$OTU %in% opp_F4_10x_F2$OTU, 1, 0)
diff_F4_10x$F3_control_10x <- ifelse(diff_F4_10x$F3_control_10x == 1 | diff_F4_10x$OTU %in% opp_F4_10x_F3$OTU, 1, 0)

# create data frame for upset
F4_control_10x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F4_control_10x, decreasing = all.do$F4_control_10x,
                                   "Unfiltered Different" = diff_F4_10x$Unfiltered_control_10x,
                                   "F1 Different" = diff_F4_10x$F1_control_10x,
                                   "F2 Different" = diff_F4_10x$F2_control_10x,
                                   "F3 Different" = diff_F4_10x$F3_control_10x)
# get rid of OTUs that are not responding in this fraction
F4_10x_check <- subset(F4_control_10x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F4_10x_check) #134
nrow(F4_10x_check[which(F4_10x_check$increasing == 1),]) #49
nrow(F4_10x_check[which(F4_10x_check$decreasing == 1),]) #85

# vector of directly affected OTUs (the same response in all fractions)
otus.F4.10x.up.direct <- F4_10x_check$OTU[which(rowSums(F4_10x_check[,-1]) == 1 & F4_10x_check$increasing == 1)] #3
otus.F4.10x.do.direct <- F4_10x_check$OTU[which(rowSums(F4_10x_check[,-1]) == 1 & F4_10x_check$decreasing == 1)] #9
otus.F4.10x.all.direct <- c(otus.F4.10x.up.direct, otus.F4.10x.do.direct) #12

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F4.10x.indirect <- F4_10x_check$OTU[which(F4_10x_check$OTU %!in% otus.F4.10x.all.direct)] #122

# count opposite
otus.F4.10x.opp <- unique(c(opp_F4_10x_U$OTU, opp_F4_10x_F1$OTU, opp_F4_10x_F2$OTU, opp_F4_10x_F3$OTU)) #23

# put into data frame to add to affected.otus
affected.F4.10x <- data.frame("Indirect" = length(otus.F4.10x.indirect),
                              "OTUs Increasing" = nrow(F4_10x_check[which(F4_10x_check$increasing == 1),]),
                              "OTUs Decreasing" = nrow(F4_10x_check[which(F4_10x_check$decreasing == 1),]),
                              "Increasing in all fractions" = length(otus.F4.10x.up.direct), 
                              "Decreasing in all fractions" = length(otus.F4.10x.do.direct),
                              "Opposite in at least one fraction" = length(otus.F4.10x.opp),
                              row.names = "F4 Control - 10x")
affected.add <- rbind(affected.add, affected.F4.10x)

# Make upset plot
F4.10x.upset <- upset(F4_10x_check, sets = c("F3.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                      keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F4.10x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F4_10x_check$OTU), "OTUs significantly affected in F4 10x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_10x_check[which(F4_10x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_10x_check[which(F4_10x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.10x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F4.10x$Increasing.in.all.fractions/affected.F4.10x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.10x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F4.10x$Decreasing.in.all.fractions/affected.F4.10x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F4.10x$Indirect, "of", length(F4_10x_check$OTU), "OTUs (", round(affected.F4.10x$Indirect/length(F4_10x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F4.10x.upset ,file="out/F4.10x.upset.pptx")

## F4 Control - 1x ----
# find OTUs that are affected in opposite direction in another fraction
opp_F4_1x = all.1x
opp_F4_1x_U <- opp_F4_1x[which(opp_F4_1x$F4_control_1x != 0 & opp_F4_1x$Unfiltered_control_1x != 0),]
opp_F4_1x_U <- opp_F4_1x_U[!(sign(opp_F4_1x_U$F4_control_1x) == sign(opp_F4_1x_U$Unfiltered_control_1x)),]

opp_F4_1x_U <- all.1x %>%
  subset(F4_control_1x != 0 & Unfiltered_control_1x != 0) %>%
  subset(!(sign(F4_control_1x) == sign(Unfiltered_control_1x)))

opp_F4_1x_F1 <- all.1x %>%
  subset(F4_control_1x != 0 & F1_control_1x != 0) %>%
  subset(!(sign(F4_control_1x) == sign(F1_control_1x))) 

opp_F4_1x_F2 <- all.1x %>%
  subset(F4_control_1x != 0 & F2_control_1x != 0) %>%
  subset(!(sign(F4_control_1x) == sign(F2_control_1x)))

opp_F4_1x_F3 <- all.1x %>%
  subset(F4_control_1x != 0 & F3_control_1x != 0) %>%
  subset(!(sign(F4_control_1x) == sign(F3_control_1x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_F4_1x = pres.nochange.1x
diff_F4_1x$Unfiltered_control_1x <- ifelse(diff_F4_1x$Unfiltered_control_1x == 1 | diff_F4_1x$OTU %in% opp_F4_1x_U$OTU, 1, 0)
diff_F4_1x$F1_control_1x <- ifelse(diff_F4_1x$F1_control_1x == 1 | diff_F4_1x$OTU %in% opp_F4_1x_F1$OTU, 1, 0)
diff_F4_1x$F2_control_1x <- ifelse(diff_F4_1x$F2_control_1x == 1 | diff_F4_1x$OTU %in% opp_F4_1x_F2$OTU, 1, 0)
diff_F4_1x$F3_control_1x <- ifelse(diff_F4_1x$F3_control_1x == 1 | diff_F4_1x$OTU %in% opp_F4_1x_F3$OTU, 1, 0)

# create data frame for upset
F4_control_1x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$F4_control_1x, decreasing = all.do$F4_control_1x,
                                  "Unfiltered Different" = diff_F4_1x$Unfiltered_control_1x,
                                  "F1 Different" = diff_F4_1x$F1_control_1x,
                                  "F2 Different" = diff_F4_1x$F2_control_1x,
                                  "F3 Different" = diff_F4_1x$F3_control_1x)
# get rid of OTUs that are not responding in this fraction
F4_1x_check <- subset(F4_control_1x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(F4_1x_check) #61
nrow(F4_1x_check[which(F4_1x_check$increasing == 1),]) #26
nrow(F4_1x_check[which(F4_1x_check$decreasing == 1),]) #35

# vector of directly affected OTUs (the same response in all fractions)
otus.F4.1x.up.direct <- F4_1x_check$OTU[which(rowSums(F4_1x_check[,-1]) == 1 & F4_1x_check$increasing == 1)] #2
otus.F4.1x.do.direct <- F4_1x_check$OTU[which(rowSums(F4_1x_check[,-1]) == 1 & F4_1x_check$decreasing == 1)] #5
otus.F4.1x.all.direct <- c(otus.F4.1x.up.direct, otus.F4.1x.do.direct) #7

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.F4.1x.indirect <- F4_1x_check$OTU[which(F4_1x_check$OTU %!in% otus.F4.1x.all.direct)] #54

# count opposite
otus.F4.1x.opp <- unique(c(opp_F4_1x_U$OTU, opp_F4_1x_F1$OTU, opp_F4_1x_F2$OTU, opp_F4_1x_F3$OTU)) #10

# put into data frame to add to affected.otus
affected.F4.1x <- data.frame("Indirect" = length(otus.F4.1x.indirect),
                             "OTUs Increasing" = nrow(F4_1x_check[which(F4_1x_check$increasing == 1),]),
                             "OTUs Decreasing" = nrow(F4_1x_check[which(F4_1x_check$decreasing == 1),]),
                             "Increasing in all fractions" = length(otus.F4.1x.up.direct), 
                             "Decreasing in all fractions" = length(otus.F4.1x.do.direct),
                             "Opposite in at least one fraction" = length(otus.F4.1x.opp),
                             row.names = "F4 Control - 1x")
affected.add <- rbind(affected.add, affected.F4.1x)

# Make upset plot
F4.1x.upset <- upset(F4_1x_check, sets = c("F3.Different", "F2.Different", "F1.Different", "Unfiltered.Different", "decreasing", "increasing"), 
                     keep.order = TRUE, order.by = "freq", text.scale = 1.5)
F4.1x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(F4_1x_check$OTU), "OTUs significantly affected in F4 1x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_1x_check[which(F4_1x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(F4_1x_check[which(F4_1x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.1x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.F4.1x$Increasing.in.all.fractions/affected.F4.1x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.F4.1x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.F4.1x$Decreasing.in.all.fractions/affected.F4.1x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.F4.1x$Indirect, "of", length(F4_1x_check$OTU), "OTUs (", round(affected.F4.1x$Indirect/length(F4_1x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= F4.1x.upset ,file="out/F4.1x.upset.pptx")

## Unfiltered Control - 100x ----
# find OTUs that are affected in opposite direction in another fraction
opp_Unfiltered_100x_F4 <- all.100x %>%
  subset(F4_control_100x != 0 & Unfiltered_control_100x != 0) %>%
  subset(!(sign(F4_control_100x) == sign(Unfiltered_control_100x)))

opp_Unfiltered_100x_F1 <- all.100x %>%
  subset(Unfiltered_control_100x != 0 & F1_control_100x != 0) %>%
  subset(!(sign(Unfiltered_control_100x) == sign(F1_control_100x))) 

opp_Unfiltered_100x_F2 <- all.100x %>%
  subset(Unfiltered_control_100x != 0 & F2_control_100x != 0) %>%
  subset(!(sign(Unfiltered_control_100x) == sign(F2_control_100x)))

opp_Unfiltered_100x_F3 <- all.100x %>%
  subset(Unfiltered_control_100x != 0 & F3_control_100x != 0) %>%
  subset(!(sign(Unfiltered_control_100x) == sign(F3_control_100x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_Unfiltered_100x = pres.nochange.100x
diff_Unfiltered_100x$F4_control_100x <- ifelse(diff_Unfiltered_100x$F4_control_100x == 1 | diff_Unfiltered_100x$OTU %in% opp_Unfiltered_100x_F4$OTU, 1, 0)
diff_Unfiltered_100x$F1_control_100x <- ifelse(diff_Unfiltered_100x$F1_control_100x == 1 | diff_Unfiltered_100x$OTU %in% opp_Unfiltered_100x_F1$OTU, 1, 0)
diff_Unfiltered_100x$F2_control_100x <- ifelse(diff_Unfiltered_100x$F2_control_100x == 1 | diff_Unfiltered_100x$OTU %in% opp_Unfiltered_100x_F2$OTU, 1, 0)
diff_Unfiltered_100x$F3_control_100x <- ifelse(diff_Unfiltered_100x$F3_control_100x == 1 | diff_Unfiltered_100x$OTU %in% opp_Unfiltered_100x_F3$OTU, 1, 0)

# create data frame for upset
Unfiltered_control_100x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$Unfiltered_control_100x, decreasing = all.do$Unfiltered_control_100x,
                                    "F4 Different" = diff_Unfiltered_100x$F4_control_100x,
                                    "F1 Different" = diff_Unfiltered_100x$F1_control_100x,
                                    "F2 Different" = diff_Unfiltered_100x$F2_control_100x,
                                    "F3 Different" = diff_Unfiltered_100x$F3_control_100x)
# get rid of OTUs that are not responding in this fraction
Unfiltered_100x_check <- subset(Unfiltered_control_100x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(Unfiltered_100x_check) #256
nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$increasing == 1),]) #130
nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$decreasing == 1),]) #126

# vector of directly affected OTUs (the same response in all fractions)
otus.Unfiltered.100x.up.direct <- Unfiltered_100x_check$OTU[which(rowSums(Unfiltered_100x_check[,-1]) == 1 & Unfiltered_100x_check$increasing == 1)] #45
otus.Unfiltered.100x.do.direct <- Unfiltered_100x_check$OTU[which(rowSums(Unfiltered_100x_check[,-1]) == 1 & Unfiltered_100x_check$decreasing == 1)] #23
otus.Unfiltered.100x.all.direct <- c(otus.Unfiltered.100x.up.direct, otus.Unfiltered.100x.do.direct) #68

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.Unfiltered.100x.indirect <- Unfiltered_100x_check$OTU[which(Unfiltered_100x_check$OTU %!in% otus.Unfiltered.100x.all.direct)] #188

# count opposite
otus.Unfiltered.100x.opp <- unique(c(opp_Unfiltered_100x_F4$OTU, opp_Unfiltered_100x_F1$OTU, opp_Unfiltered_100x_F2$OTU, opp_Unfiltered_100x_F3$OTU)) #29

# put into data frame to add to affected.otus
affected.Unfiltered.100x <- data.frame("Indirect" = length(otus.Unfiltered.100x.indirect),
                               "OTUs Increasing" = nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$increasing == 1),]),
                               "OTUs Decreasing" = nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$decreasing == 1),]),
                               "Increasing in all fractions" = length(otus.Unfiltered.100x.up.direct), 
                               "Decreasing in all fractions" = length(otus.Unfiltered.100x.do.direct),
                               "Opposite in at least one fraction" = length(otus.Unfiltered.100x.opp),
                               row.names = "Unfiltered Control - 100x")
affected.add <- rbind(affected.add, affected.Unfiltered.100x)

# Make upset plot
Unfiltered.100x.upset <- upset(Unfiltered_100x_check, sets = c("F4.Different","F3.Different", "F2.Different", "F1.Different", "decreasing", "increasing"), 
                       keep.order = TRUE, order.by = "freq", text.scale = 1.5)
Unfiltered.100x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(Unfiltered_100x_check$OTU), "OTUs significantly affected in Unfiltered 100x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_100x_check[which(Unfiltered_100x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.100x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.Unfiltered.100x$Increasing.in.all.fractions/affected.Unfiltered.100x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.100x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.Unfiltered.100x$Decreasing.in.all.fractions/affected.Unfiltered.100x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.Unfiltered.100x$Indirect, "of", length(Unfiltered_100x_check$OTU), "OTUs (", round(affected.Unfiltered.100x$Indirect/length(Unfiltered_100x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= Unfiltered.100x.upset ,file="out/Unfiltered.100x.upset.pptx")

## Unfiltered Control - 10x ----
# find OTUs that are affected in opposite direction in another fraction
opp_Unfiltered_10x_F4 <- all.10x %>%
  subset(F4_control_10x != 0 & Unfiltered_control_10x != 0) %>%
  subset(!(sign(F4_control_10x) == sign(Unfiltered_control_10x)))

opp_Unfiltered_10x_F1 <- all.10x %>%
  subset(Unfiltered_control_10x != 0 & F1_control_10x != 0) %>%
  subset(!(sign(Unfiltered_control_10x) == sign(F1_control_10x))) 

opp_Unfiltered_10x_F2 <- all.10x %>%
  subset(Unfiltered_control_10x != 0 & F2_control_10x != 0) %>%
  subset(!(sign(Unfiltered_control_10x) == sign(F2_control_10x)))

opp_Unfiltered_10x_F3 <- all.10x %>%
  subset(Unfiltered_control_10x != 0 & F3_control_10x != 0) %>%
  subset(!(sign(Unfiltered_control_10x) == sign(F3_control_10x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_Unfiltered_10x = pres.nochange.10x
diff_Unfiltered_10x$F4_control_10x <- ifelse(diff_Unfiltered_10x$F4_control_10x == 1 | diff_Unfiltered_10x$OTU %in% opp_Unfiltered_10x_F4$OTU, 1, 0)
diff_Unfiltered_10x$F1_control_10x <- ifelse(diff_Unfiltered_10x$F1_control_10x == 1 | diff_Unfiltered_10x$OTU %in% opp_Unfiltered_10x_F1$OTU, 1, 0)
diff_Unfiltered_10x$F2_control_10x <- ifelse(diff_Unfiltered_10x$F2_control_10x == 1 | diff_Unfiltered_10x$OTU %in% opp_Unfiltered_10x_F2$OTU, 1, 0)
diff_Unfiltered_10x$F3_control_10x <- ifelse(diff_Unfiltered_10x$F3_control_10x == 1 | diff_Unfiltered_10x$OTU %in% opp_Unfiltered_10x_F3$OTU, 1, 0)

# create data frame for upset
Unfiltered_control_10x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$Unfiltered_control_10x, decreasing = all.do$Unfiltered_control_10x,
                                            "F4 Different" = diff_Unfiltered_10x$F4_control_10x,
                                            "F1 Different" = diff_Unfiltered_10x$F1_control_10x,
                                            "F2 Different" = diff_Unfiltered_10x$F2_control_10x,
                                            "F3 Different" = diff_Unfiltered_10x$F3_control_10x)
# get rid of OTUs that are not responding in this fraction
Unfiltered_10x_check <- subset(Unfiltered_control_10x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(Unfiltered_10x_check) #175
nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$increasing == 1),]) #102
nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$decreasing == 1),]) #73

# vector of directly affected OTUs (the same response in all fractions)
otus.Unfiltered.10x.up.direct <- Unfiltered_10x_check$OTU[which(rowSums(Unfiltered_10x_check[,-1]) == 1 & Unfiltered_10x_check$increasing == 1)] #10
otus.Unfiltered.10x.do.direct <- Unfiltered_10x_check$OTU[which(rowSums(Unfiltered_10x_check[,-1]) == 1 & Unfiltered_10x_check$decreasing == 1)] #13
otus.Unfiltered.10x.all.direct <- c(otus.Unfiltered.10x.up.direct, otus.Unfiltered.10x.do.direct) #23

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.Unfiltered.10x.indirect <- Unfiltered_10x_check$OTU[which(Unfiltered_10x_check$OTU %!in% otus.Unfiltered.10x.all.direct)] #152

# count opposite
otus.Unfiltered.10x.opp <- unique(c(opp_Unfiltered_10x_F4$OTU, opp_Unfiltered_10x_F1$OTU, opp_Unfiltered_10x_F2$OTU, opp_Unfiltered_10x_F3$OTU)) #31

# put into data frame to add to affected.otus
affected.Unfiltered.10x <- data.frame("Indirect" = length(otus.Unfiltered.10x.indirect),
                                       "OTUs Increasing" = nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$increasing == 1),]),
                                       "OTUs Decreasing" = nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$decreasing == 1),]),
                                       "Increasing in all fractions" = length(otus.Unfiltered.10x.up.direct), 
                                       "Decreasing in all fractions" = length(otus.Unfiltered.10x.do.direct),
                                       "Opposite in at least one fraction" = length(otus.Unfiltered.10x.opp),
                                       row.names = "Unfiltered Control - 10x")
affected.add <- rbind(affected.add, affected.Unfiltered.10x)

# Make upset plot
Unfiltered.10x.upset <- upset(Unfiltered_10x_check, sets = c("F4.Different","F3.Different", "F2.Different", "F1.Different", "decreasing", "increasing"), 
                               keep.order = TRUE, order.by = "freq", text.scale = 1.5)
Unfiltered.10x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(Unfiltered_10x_check$OTU), "OTUs significantly affected in Unfiltered 10x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_10x_check[which(Unfiltered_10x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.10x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.Unfiltered.10x$Increasing.in.all.fractions/affected.Unfiltered.10x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.10x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.Unfiltered.10x$Decreasing.in.all.fractions/affected.Unfiltered.10x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.Unfiltered.10x$Indirect, "of", length(Unfiltered_10x_check$OTU), "OTUs (", round(affected.Unfiltered.10x$Indirect/length(Unfiltered_10x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= Unfiltered.10x.upset ,file="out/Unfiltered.10x.upset.pptx")

## Unfiltered Control - 1x ----
# find OTUs that are affected in opposite direction in another fraction
opp_Unfiltered_1x_F4 <- all.1x %>%
  subset(F4_control_1x != 0 & Unfiltered_control_1x != 0) %>%
  subset(!(sign(F4_control_1x) == sign(Unfiltered_control_1x)))

opp_Unfiltered_1x_F1 <- all.1x %>%
  subset(Unfiltered_control_1x != 0 & F1_control_1x != 0) %>%
  subset(!(sign(Unfiltered_control_1x) == sign(F1_control_1x))) 

opp_Unfiltered_1x_F2 <- all.1x %>%
  subset(Unfiltered_control_1x != 0 & F2_control_1x != 0) %>%
  subset(!(sign(Unfiltered_control_1x) == sign(F2_control_1x)))

opp_Unfiltered_1x_F3 <- all.1x %>%
  subset(Unfiltered_control_1x != 0 & F3_control_1x != 0) %>%
  subset(!(sign(Unfiltered_control_1x) == sign(F3_control_1x)))

# add response to no.change so that an OTU is "1" if it is not responding or responding differently in another fraction
diff_Unfiltered_1x = pres.nochange.1x
diff_Unfiltered_1x$F4_control_1x <- ifelse(diff_Unfiltered_1x$F4_control_1x == 1 | diff_Unfiltered_1x$OTU %in% opp_Unfiltered_1x_F4$OTU, 1, 0)
diff_Unfiltered_1x$F1_control_1x <- ifelse(diff_Unfiltered_1x$F1_control_1x == 1 | diff_Unfiltered_1x$OTU %in% opp_Unfiltered_1x_F1$OTU, 1, 0)
diff_Unfiltered_1x$F2_control_1x <- ifelse(diff_Unfiltered_1x$F2_control_1x == 1 | diff_Unfiltered_1x$OTU %in% opp_Unfiltered_1x_F2$OTU, 1, 0)
diff_Unfiltered_1x$F3_control_1x <- ifelse(diff_Unfiltered_1x$F3_control_1x == 1 | diff_Unfiltered_1x$OTU %in% opp_Unfiltered_1x_F3$OTU, 1, 0)

# create data frame for upset
Unfiltered_control_1x_check <- data.frame(OTU = all.df$OTU, increasing = all.up$Unfiltered_control_1x, decreasing = all.do$Unfiltered_control_1x,
                                            "F4 Different" = diff_Unfiltered_1x$F4_control_1x,
                                            "F1 Different" = diff_Unfiltered_1x$F1_control_1x,
                                            "F2 Different" = diff_Unfiltered_1x$F2_control_1x,
                                            "F3 Different" = diff_Unfiltered_1x$F3_control_1x)
# get rid of OTUs that are not responding in this fraction
Unfiltered_1x_check <- subset(Unfiltered_control_1x_check, increasing + decreasing != 0)


# Count total OTUs responding in this comparison
nrow(Unfiltered_1x_check) #72
nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$increasing == 1),]) #33
nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$decreasing == 1),]) #39

# vector of directly affected OTUs (the same response in all fractions)
otus.Unfiltered.1x.up.direct <- Unfiltered_1x_check$OTU[which(rowSums(Unfiltered_1x_check[,-1]) == 1 & Unfiltered_1x_check$increasing == 1)] #2
otus.Unfiltered.1x.do.direct <- Unfiltered_1x_check$OTU[which(rowSums(Unfiltered_1x_check[,-1]) == 1 & Unfiltered_1x_check$decreasing == 1)] #9
otus.Unfiltered.1x.all.direct <- c(otus.Unfiltered.1x.up.direct, otus.Unfiltered.1x.do.direct) #11

#vector of indirectly affected otus
`%!in%` <- Negate(`%in%`)
otus.Unfiltered.1x.indirect <- Unfiltered_1x_check$OTU[which(Unfiltered_1x_check$OTU %!in% otus.Unfiltered.1x.all.direct)] #61

# count opposite
otus.Unfiltered.1x.opp <- unique(c(opp_Unfiltered_1x_F4$OTU, opp_Unfiltered_1x_F1$OTU, opp_Unfiltered_1x_F2$OTU, opp_Unfiltered_1x_F3$OTU)) #14

# put into data frame to add to affected.otus
affected.Unfiltered.1x <- data.frame("Indirect" = length(otus.Unfiltered.1x.indirect),
                                       "OTUs Increasing" = nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$increasing == 1),]),
                                       "OTUs Decreasing" = nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$decreasing == 1),]),
                                       "Increasing in all fractions" = length(otus.Unfiltered.1x.up.direct), 
                                       "Decreasing in all fractions" = length(otus.Unfiltered.1x.do.direct),
                                       "Opposite in at least one fraction" = length(otus.Unfiltered.1x.opp),
                                       row.names = "Unfiltered Control - 1x")
affected.add <- rbind(affected.add, affected.Unfiltered.1x)

# Make upset plot
Unfiltered.1x.upset <- upset(Unfiltered_1x_check, sets = c("F4.Different","F3.Different", "F2.Different", "F1.Different", "decreasing", "increasing"), 
                               keep.order = TRUE, order.by = "freq", text.scale = 1.5)
Unfiltered.1x.upset
# UpSetR does not have functionality for adding labels. Add grid.text to avoid confusion when looking at upsets
grid.text(paste(length(Unfiltered_1x_check$OTU), "OTUs significantly affected in Unfiltered 1x compared to control"), x = 0.5, y=0.96, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$increasing == 1),]), "OTUs increasing"), x = 0.5, y=0.93, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(nrow(Unfiltered_1x_check[which(Unfiltered_1x_check$decreasing == 1),]), "OTUs decreasing"), x = 0.5, y=0.90, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.1x.up.direct), "OTUs increasing in all fractions (", 
                round((affected.Unfiltered.1x$Increasing.in.all.fractions/affected.Unfiltered.1x$OTUs.Increasing)*100, digits = 2), "% of increasing)"), x = 0.5, y=0.87, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(length(otus.Unfiltered.1x.do.direct), "OTUs decreasing in all fractions (", 
                round((affected.Unfiltered.1x$Decreasing.in.all.fractions/affected.Unfiltered.1x$OTUs.Decreasing)*100, digits = 2), "% of decreasing)"), x = 0.5, y=0.84, gp=gpar(fontsize=16), just = "left") 
grid.text(paste(affected.Unfiltered.1x$Indirect, "of", length(Unfiltered_1x_check$OTU), "OTUs (", round(affected.Unfiltered.1x$Indirect/length(Unfiltered_1x_check$OTU)*100, digits = 2), 
                "%) due to indirect effects"), x = 0.5, y=0.81, gp=gpar(fontsize=16), just = "left") 

#also save as .pptx for better ability to manipulate labels, etc
graph2ppt(x= Unfiltered.1x.upset ,file="out/Unfiltered.1x.upset.pptx")


