## ----setup, warning = FALSE, message=FALSE------------------------------------
library(DiSC)
require(GUniFrac)

## ----input--------------------------------------------------------------------
data("sim_data", package = "DiSC")

count_matrix <- sim_data$count_matrix
meta_cell <- sim_data$meta_cell
gene_index <- sim_data$gene_index
meta_ind <- sim_data$meta_ind

dim(count_matrix)
dim(meta_cell)
dim(meta_ind)
count_matrix[1:3,1:6]
head(meta_cell)
head(meta_ind)

## ----rarefaction--------------------------------------------------------------
read_depth = rowSums(count_matrix)
discarding_flag = FALSE
if(discarding_flag){
  # when the sequencing depths of some samples are shallow
  depth = min(5000, round(quantile(read_depth, probs = 0.1)))
  cat("read depth filtering with threshold", depth, "\n")
  count_matrix_temp <- count_matrix[, which(read_depth >= depth)]
  meta_cell_temp <- meta_cell %>% 
    filter(cell_id %in% colnames(count_matrix_temp))
  rarefy_mat <- t(GUniFrac::Rarefy(t(count_matrix_temp))$otu.tab.rff)
} else {
  # when all samples are sufficiently sequenced
  rarefy_mat <- t(GUniFrac::Rarefy(t(count_matrix))$otu.tab.rff)
}


## -----------------------------------------------------------------------------
obj1 <- DiSC(ct.mat = rarefy_mat, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 999,
             features = c('prev', 'nzm', 'nzsd'), verbose = TRUE)
# Type I error
mean(obj1$p.raw[gene_index$EE_index] <= 0.05)
# Power
mean(obj1$p.raw[gene_index$mean_index] <= 0.05)
mean(obj1$p.raw[gene_index$var_index] <= 0.05)
mean(obj1$p.raw[gene_index$mean_var_index] <= 0.05)
# False discovery proportion
sum(obj1$p.adj.fdr[gene_index$EE_index] <= 0.10)/
  sum(obj1$p.adj.fdr <= 0.10)
# Number of positive discoveries
sum(obj1$p.adj.fdr[gene_index$mean_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$var_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$mean_var_index] <= 0.10)

t <- proc.time()
obj2 <- DiSC(ct.mat = rarefy_mat, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 99,
             features = c('prev', 'nzm', 'nzsd'), verbose = FALSE)
print(t - proc.time())
obj1$p.adj.fdr[123]
obj2$p.adj.fdr[123]
sum(obj1$p.adj.fdr[gene_index$mean_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$var_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$mean_var_index] <= 0.10)

