# Introduction ####
# This script is used to compute gene coexpression and Gene Coexpression Contribution (GCC) scores 
# (Fulcher & Fornito, PNAS 2016). We compute mean gene coexpression for each of 19 brain regions of interest, 
# and correct for spatial correlations in the data. GCC scores are then computed in order to identify genes that 
# significantly contribute to higher gene coexpression observed amongst pairs of brain regions that show increased
# weighted degree following opioid dependence (change from naive to 24 h state) compared to pairs of brain regions that
# show decreased weighted degree following opioid dependence.

# 1. Read in and prepare the datasets for analysis ####
# Prepare the gene expression data
# read in the Allen gene expression dataset
allGE <- read.csv("~/Downloads/AllGeneExpression_NormalizedFirst.csv", row.names=1,)
# select the data frame rows corresponding to brain regions of interest
GE <- as.data.frame(allGE[c(2,3,7,9,31,36,4,164,24,22,28,41,92,74,146,117,167,168,212), ])
# remove genes with missing values (NA)
GE <- GE[ , colSums(is.na(GE)) == 0]
# transpose the dataframe so brain regions are columns
tGE <- as.data.frame(t(GE))

# Prepare in the distance data
# read in the Allen distance dataset
all_dist <- read.csv("~/Downloads/nature13186-s5.csv", row.names=1,)
# select the data frame rows corresponding to ipsilateral data for the brain regions of interest
dist.mat <- as.matrix(all_dist[c("ACAd_ipsi","ACAv_ipsi", "AId_ipsi", "AIv_ipsi", "CLA_ipsi", "CP_ipsi", "ACB_ipsi", "SI_ipsi", "BST_ipsi", "BLA_ipsi", "CEA_ipsi", "DG_ipsi", "MH_ipsi", "LH_ipsi","PVT_ipsi", "PAG_ipsi", "SNc_ipsi", "SNr_ipsi", "VTA_ipsi"), 
                               c("ACAd_ipsi","ACAv_ipsi", "AId_ipsi", "AIv_ipsi", "CLA_ipsi", "CP_ipsi", "ACB_ipsi", "SI_ipsi", "BST_ipsi", "BLA_ipsi", "CEA_ipsi", "DG_ipsi", "MH_ipsi", "LH_ipsi","PVT_ipsi", "PAG_ipsi", "SNc_ipsi", "SNr_ipsi", "VTA_ipsi")])
# convert the distances to mm
dist.mat <- dist.mat / 1000
# preserve the lower triangle of the matrix
dist <- dist.mat[lower.tri(dist.mat)]

# Prepare the FOS data
# read in data for naive and 24 h conditions
Fos_naive <- read.csv("~/Desktop/FCbyRegion_Naive.csv")
Fos_24h <- read.csv("~/Desktop/FCbyRegion_24h.csv")
# obtain matrix of pairwise correlations in FOS data
corFos_naive<-cor(Fos_naive)
corFos_24h <- cor(Fos_24h)
# replace negative correlations with 0 to obtain edge weights
corFos_naive[corFos_naive<0] <- 0
corFos_24h[corFos_24h<0] <- 0
# preserve the lower triangle of the matrix
corFos_naive <- corFos_naive[lower.tri(corFos_naive)]
corFos_24h <- corFos_24h[lower.tri(corFos_24h)]

# compute the change in edge weight from naive to 24 h
change_weight <- corFos_24h - corFos_naive

# 2. Compute gene coexpression across all genes, between every pair of brain regions ####
# create a matrix of correlated gene expression across all genes, between every pair of brain regions
cor_tGE.mat <- cor(tGE)
# preserve the lower triangle of the matrix
cor_tGE <- cor_tGE.mat[lower.tri(cor_tGE.mat)]
# bind gene expression and distance data
DistvGEcor <- as.data.frame(cbind(dist, cor_tGE))
# graph distance vs gene coexpression
plot(dist, cor_tGE)
# I wrote the distance and gene expression out to a file, and then used Prism software
# to fit an exponential decay curve to the data
write.csv(DistvGEcor, "~/Desktop/DistvGEcor.csv")

# 3. Compute spatially-corrected mean gene coexpression per brain region ####
# compute the spatially-corrected gene coexpression for every pair of brain regions using the exponential decay equation
# and distance data
n=171
DistvGEcor$cor_tGE_corrected <- as.numeric(lapply(1:n, function(x) (DistvGEcor[x,2] - (1.2201*(2.71828^(-0.6008*(DistvGEcor[x,1]))) - 0.1491))))
# visualize the effect of the spatial correction on the relationship between distance and gene coexpression
plot(DistvGEcor$dist, DistvGEcor$cor_tGE_corrected)

# Compute uncorrected mean gene coexpression
cor_tGE.mat.nd <- cor_tGE.mat
# remove the diagonal (self-correlations)
diag(cor_tGE.mat.nd)=NA
cor_tGE.mat.nd <- matrix(t(cor_tGE.mat.nd)[which(!is.na(cor_tGE.mat.nd))],nrow=18,ncol=19)
# get mean gene coexpression per brain region
meanGE <- apply(cor_tGE.mat.nd, 2, mean)

# Compute spatially-corrected mean gene coexpression
# perform spatial correction over the matrix of gene coexpression values
correctedGE <- sapply(1:361, function(x) (cor_tGE.mat[[x]] - (1.2201*(2.71828^(-0.6008*(dist.mat[[x]]))) - 0.1491)))
# convert the output into a matrix
correctedGE.mat <- matrix(correctedGE, 19, 19)
# remove the diagonal (self-correlations)
diag(correctedGE.mat)=NA
correctedGE.mat <- matrix(t(correctedGE.mat)[which(!is.na(correctedGE.mat))],nrow=18,ncol=19)
# get mean gene coexpression per brain region
meanGEcorrect <- apply(correctedGE.mat, 2, mean)

# 4. Compute Gene Coexpression Contribution (GCC) scores ####
# compute the spatial correction factor for each pair of brain regions
n=171
spatial_correct_factor <- lapply(1:n, function(x) (1.2201*(2.71828^(-0.6008*(DistvGEcor[x,1]))) - 0.1491))
spatial_correct_factor <- as.data.frame(spatial_correct_factor)

# z-score the genes within each brain region
# z-score relative to all genes in each brain region
zscore_byregion<-apply(tGE, 2, scale)
# convert the output to a dataframe; preserve the gene names
zscore_byregion<-as.data.frame(zscore_byregion, row.names = row.names(tGE))

# For every gene, take the product of z-scored expression for every pair of brain regions
# set up a function for all pairwise combinations of brain regions
n=ncol(zscore_byregion)
combb=combn(n,2)
combb=cbind(combb, sapply(1:n, function(i) rep(i,2)))
# take the product of gene a in every pair of regions i and j
zscore_prods=apply(zscore_byregion, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
# convert to dataframe
zscore_prods<-as.data.frame(t(zscore_prods))
# exclude the self-products
zscore_prods <- zscore_prods[-c(172:190)]

# Subtract the spatially-corrected correlation value for every pair of brain regions to get the gene coexpression contribution!
n=171
GCC<-lapply(1:n, function(x) (zscore_prods[,x] - spatial_correct_factor[,x]))
GCC<-as.data.frame(t(GCC))
colnames(GCC) <- row.names(tGE)

# subset gene expression data into increasing and decreasing edge weight groups
inc_weightGCC <- GCC[c(7,12,23:24,26,29,31,38:40,42:45,47,53,57,60:64,67:69,71:74,76,82,84:85,89,94,96:97,99,101,106:108,110,112,122,127:128,131,136,139,166,168),]
dec_weightGCC <- GCC[-c(7,12,23:24,26,29,31,38:40,42:45,47,53,57,60:64,67:69,71:74,76,82,84:85,89,94,96:97,99,101,106:108,110,112,122,127:128,131,136,139,166,168),]

# perform one-sided t-test for each gene
n=19616
ttest_GCC <- lapply(1:n, function(x) t.test(inc_weightGCC[,x], dec_weightGCC[,x], alternative="greater"))
ttest_GCC<-setNames(ttest_GCC,row.names(tGE))

# compile the t-test results for each gene into a dataframe
GCC_est<-lapply(ttest_GCC, function(x) x$estimate)
GCC_est<-as.data.frame(GCC_est)
GCC_est<-as.data.frame(t(GCC_est))
GCC_est$Gene<-colnames(GE)

GCC_results<-lapply(ttest_GCC, function(x) x$statistic)
GCC_results<-as.data.frame(GCC_results)
GCC_results<-as.data.frame(t(GCC_results))
GCC_results$Gene<-colnames(GE)

GCC_pval<-lapply(ttest_GCC, function(x) x$p.value)
GCC_pval<-as.data.frame(GCC_pval)
GCC_pval<-as.data.frame(t(GCC_pval))
GCC_pval$Gene<-colnames(GE)

# adjust p-values using FDR
cGCC_pval<-p.adjust(GCC_pval$V1, method = "fdr", n = length(GCC_pval$V1))
cGCC_pval<-as.data.frame(cGCC_pval)
cGCC_pval$Gene<-colnames(GE)
totalGCC_results <- merge(GCC_est, GCC_results,by="Gene")
totalGCC_results <- merge(totalGCC_results, GCC_pval,by="Gene")
totalGCC_results <- merge(totalGCC_results, cGCC_pval,by="Gene")
totalGCC_results <- totalGCC_results[order(totalGCC_results$cGCC_pval),]

# write out the results to a single file
write.csv(totalGCC_results, "~Desktop/totalGCC_results.csv")