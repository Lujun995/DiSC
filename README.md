# DiSC: A statistical tool for differential expression analyis of individual level single-cell RNA-Seq data

## Installation 

 
```R
devtools::install_github("Lujun995/DiSC")
library(DiSC)
```

## Vignette

```R
vignette("DiSC")
```


## Usage


```R
data("sim_data", package = "DiSC")

count_matrix <- sim_data$count_matrix
meta_cell <- sim_data$meta_cell
gene_index <- sim_data$gene_index
meta_ind <- sim_data$meta_ind

set.seed(seed = 123456)
t <- proc.time()
obj1 <- DiSC(data.mat = count_matrix, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 999,
             features = c('prev', 'nzm', 'nzsd'), verbose = TRUE,
             sequencing.data = TRUE)
print("Computational time for DiSC:")
print(proc.time() - t)
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
```


## Citation

Zhang, L., Yang, L., Ren, Y., Zhang, S., Guan, W., & Chen, J. (Submitted to Bioinformatics): DiSC: a Statistical Tool for Fast Differential Expression Analysis of Individual-level Single-cell RNA-seq Data.
