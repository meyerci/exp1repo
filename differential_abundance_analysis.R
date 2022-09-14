library(phyloseq)
library(readr)
library(dplyr)
library(phytools)
library(emmeans)
library(tibble)
library(lme4)
library(reshape2)

rm(list = names(.GlobalEnv)[grep("tmp",names(.GlobalEnv))])

# select data from the chosen barcode
tmp_ps0 = ps_16S_most_abund.1

# subset data by fraction
tmp_ps_U = subset_samples(tmp_ps0,tmp_ps0@sam_data$fraction  == "Unfiltered")
tmp_ps_F1 = subset_samples(tmp_ps0,tmp_ps0@sam_data$fraction  == "F1")
tmp_ps_F2 = subset_samples(tmp_ps0,tmp_ps0@sam_data$fraction  == "F2")
tmp_ps_F3 = subset_samples(tmp_ps0,tmp_ps0@sam_data$fraction  == "F3")
tmp_ps_F4 = subset_samples(tmp_ps0,tmp_ps0@sam_data$fraction  == "F4")

# select which fraction to analyze 
##### ALSO CHANGE THE NAME OF THE DATA FRAME AND MODEL AT END OF LOOP
tmp_ps = tmp_ps_F4

# treatments
a = tibble("id"= tmp_ps@sam_data$id,
           "dose"= as.character(tmp_ps@sam_data$dose))
a = factor(a$dose, levels = c("Control", "1x", "10x", "100x"))
# offset
o = log(sample_sums(tmp_ps))
# random effect
z <- as.factor(tmp_ps@sam_data$id)


# r multiple comparaison loop with model

glmPLN.sum.global = data.frame()
glmPLN.pairwise.global = data.frame()


for (i in 1:ntaxa(tmp_ps)) {
    # select one OTU
  OTU = taxa_names(tmp_ps)[i]
    # response variable
  y = as.vector(tmp_ps@otu_table[OTU,]@.Data)
    
  tryCatch({
      ### model
    glmPLN <- glmer(y ~ -1 + a + (1 | z),
                     family='poisson', offset = o)
      
    glmPLN.sum = summary(glmPLN)$coefficients
    glmPLN.sum = tibble("OTU"= OTU,
                          "treatment"=rownames(glmPLN.sum),
                          as_tibble(glmPLN.sum))
    glmPLN.sum
    glmPLN.sum.global = rbind(glmPLN.sum.global,glmPLN.sum)
      ### multiple comparaison
    glmPLN.pairwise = emmeans(glmPLN,pairwise~a,adjust="tukey")
      # select p value
    glmPLN.pairwise.sum = summary(glmPLN.pairwise)
    glmPLN.pairwise.sum = glmPLN.pairwise.sum[["contrasts"]]
      # extract summary
    tmp_df = glmPLN.pairwise.sum
      # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    tmp_df = tmp_df[tmp_df$a == "Control" | tmp_df$b == "Control",]
    tmp_df = cbind("OTU"=OTU,tmp_df)
      # extract results in data frame
    glmPLN.pairwise.global = rbind(glmPLN.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmPLN,glmPLN.sum)
}
glmPLN.model.global.F4 = glmPLN.sum.global
glmPLN.pairwise.global.F4 = glmPLN.pairwise.global

# adjust p values for all OTUs
glmPLN.pairwise.global.U[,"p.adjust"] <- p.adjust(glmPLN.pairwise.global.U$p.value,"fdr")
glmPLN.pairwise.global.F1[,"p.adjust"] <- p.adjust(glmPLN.pairwise.global.F1$p.value,"fdr")
glmPLN.pairwise.global.F2[,"p.adjust"] <- p.adjust(glmPLN.pairwise.global.F2$p.value,"fdr")
glmPLN.pairwise.global.F3[,"p.adjust"] <- p.adjust(glmPLN.pairwise.global.F3$p.value,"fdr")
glmPLN.pairwise.global.F4[,"p.adjust"] <- p.adjust(glmPLN.pairwise.global.F4$p.value,"fdr")

# Filter significant OTUs ----
## U ----
## number of significant p val <= 0.05 before and after adjustment
table(glmPLN.pairwise.global.U$p.value <= 0.05) #577
table(glmPLN.pairwise.global.U$p.adjust <= 0.05) #503

## number of OTUs with a pval <= 0.05 before and after adjustment
otu.U <- unique(glmPLN.pairwise.global.U$OTU[glmPLN.pairwise.global.U$p.value <= 0.05])
length(otu.U) #384 OTUs significant without p value adjustment
adj.otu.U <- unique(glmPLN.pairwise.global.U$OTU[glmPLN.pairwise.global.U$p.adjust <= 0.05])
length(adj.otu.U) #346 OTUs significat with p value adjustment
glmPLN.pairwise.global.U.sig <- glmPLN.pairwise.global.U[glmPLN.pairwise.global.U$p.adjust <=0.05,]

## F1 ----
## number of significant p val <= 0.05 before and after adjustment
table(glmPLN.pairwise.global.F1$p.value <= 0.05) #565
table(glmPLN.pairwise.global.F1$p.adjust <= 0.05) #490

## number of OTUs with a pval <= 0.05 before and after adjustment
otu.F1 <- unique(glmPLN.pairwise.global.F1$OTU[glmPLN.pairwise.global.F1$p.value <= 0.05])
length(otu.F1) #364 OTUs significant without p value adjustment
adj.otu.F1 <- unique(glmPLN.pairwise.global.F1$OTU[glmPLN.pairwise.global.F1$p.adjust <= 0.05])
length(adj.otu.F1) #329 OTUs significat with p value adjustment
glmPLN.pairwise.global.F1.sig <- glmPLN.pairwise.global.F1[glmPLN.pairwise.global.F1$p.adjust <=0.05,]

## F2 ----
## number of significant p val <= 0.05 before and after adjustment
table(glmPLN.pairwise.global.F2$p.value <= 0.05) #353
table(glmPLN.pairwise.global.F2$p.adjust <= 0.05) #280

## number of OTUs with a pval <= 0.05 before and after adjustment
otu.F2 <- unique(glmPLN.pairwise.global.F2$OTU[glmPLN.pairwise.global.F2$p.value <= 0.05])
length(otu.F2) #234 OTUs significant without p value adjustment
adj.otu.F2 <- unique(glmPLN.pairwise.global.F2$OTU[glmPLN.pairwise.global.F2$p.adjust <= 0.05])
length(adj.otu.F2) #195 OTUs significat with p value adjustment
glmPLN.pairwise.global.F2.sig <- glmPLN.pairwise.global.F2[glmPLN.pairwise.global.F2$p.adjust <=0.05,]

## F3 ----
## number of significant p val <= 0.05 before and after adjustment
table(glmPLN.pairwise.global.F3$p.value <= 0.05) #564
table(glmPLN.pairwise.global.F3$p.adjust <= 0.05) #476

## number of OTUs with a pval <= 0.05 before and after adjustment
otu.F3 <- unique(glmPLN.pairwise.global.F3$OTU[glmPLN.pairwise.global.F3$p.value <= 0.05])
length(otu.F3) #352 OTUs significant without p value adjustment
adj.otu.F3 <- unique(glmPLN.pairwise.global.F3$OTU[glmPLN.pairwise.global.F3$p.adjust <= 0.05])
length(adj.otu.F3) #317 OTUs significat with p value adjustment
glmPLN.pairwise.global.F3.sig <- glmPLN.pairwise.global.F3[glmPLN.pairwise.global.F3$p.adjust <=0.05,]

## F4 ----
## number of significant p val <= 0.05 before and after adjustment
table(glmPLN.pairwise.global.F4$p.value <= 0.05) #395
table(glmPLN.pairwise.global.F4$p.adjust <= 0.05) #338

## number of OTUs with a pval <= 0.05 before and after adjustment
otu.F4 <- unique(glmPLN.pairwise.global.F4$OTU[glmPLN.pairwise.global.F4$p.value <= 0.05])
length(otu.F4) #247 OTUs significant without p value adjustment
adj.otu.F4 <- unique(glmPLN.pairwise.global.F4$OTU[glmPLN.pairwise.global.F4$p.adjust <= 0.05])
length(adj.otu.F4) #211 OTUs significat with p value adjustment
glmPLN.pairwise.global.F4.sig <- glmPLN.pairwise.global.F4[glmPLN.pairwise.global.F4$p.adjust <=0.05,]

# cast pvalues
contrasts.U <- glmPLN.pairwise.global.U[,c(10,1,2)]
# numeric variable needs to be named "value" 
colnames(contrasts.U) <- c("value", "OTU_names", "contrast")
#contrasts.glm.CBFP.T2 <- subset(contrasts.glm.CBFP.T2, (contrasts.glm.CBFP.T2$OTU_names %in% BFPOTUs.T2net.sig))
head(contrasts.U)
str(contrasts.U)
contrasts.U <- dcast(contrasts.U, contrast ~ OTU_names, value="value")
str(contrasts.U)
rownames(contrasts.U) <- contrasts.U$contrast
contrasts.U$contrast <- NULL

# keep OTUs with at least one contrast <0.05 
contrasts.U.sub <- contrasts.U[,colSums(contrasts.U <0.05) >= 1]
dim(contrasts.U.sub)
head(contrasts.U.sub)
str(contrasts.U.sub)

contrasts.otu.U <- data.frame(t(contrasts.U.sub))
colnames(contrasts.otu.U)
str(contrasts.otu.U)
bin.sig.otu.U = contrasts.otu.U
bin.sig.otu.U[bin.sig.otu.U <0.05] <- 1
bin.sig.otu.U[bin.sig.otu.U !=1] <- 0
head(bin.sig.otu.U)
colSums(bin.sig.otu.U)
#Control...100x  Control...10x   Control...1x 
#330             261             256 

sig.otu.U.new <- rownames(bin.sig.otu.U) #346

tmp_ps_U 

relabund_T2 = as_tibble(tmp1T2 %>% group_by(OTUs,mesh_size_um, group, gp_sum) %>% summarise(sum=sum(Abundance), 
                                                                                            avg=mean(Abundance)) %>% summarise(rel_abund=(sum*100/gp_sum), avg=avg)) 

