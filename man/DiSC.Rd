\name{DiSC}
\alias{DiSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DiSC: A statistical tool for differential expression analyis of individual level single-cell RNA-Seq data
}
\description{
A statistical tool for differential expression analyis of individual level single-cell RNA-Seq data
}
\usage{
DiSC(data.mat, cell.ind, metadata, outcome, covariates = NULL,
     cell.id = "cell_id",
     individual.id = "individual", perm.no = 999,
     features = c('prev', 'nzm', 'nzsd'), verbose = TRUE,
     sequencing.data = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.mat}{
**Count** matrix for single cell RNA sequencing data, or other single cell data matrices. Genes / variables are in rows, **column names are cell ids** and row names are gene / variable names.
}
  \item{cell.ind}{
A data frame includes two columns of cell ids and individual ids. It links cell ids to individual ids.
}
  \item{metadata}{
A data frame includes a column of individual id, a column of an outcome of interest and columns of covariates if applicable.
}
  \item{outcome}{
A character string which is the variable name of outcome in \code{metadata}.
}
  \item{covariates}{
A character string or vector of character strings which are the variable names of covariates in \code{metadata}. Default: \code{NULL}.
}
  \item{cell.id}{
A character string which is the variable name of cell ids in \code{cell.ind}.}
  \item{individual.id}{
A character string which is the variable name of the individual id variables in \code{cell.ind} and \code{metadata}.}
  \item{perm.no}{
An integer, number of permutations used. Default: \code{999}. It can be reduced to \code{99} if adjusted P-values (false discovery rate) are the only interest.}
  \item{features}{
Vector of a character string. Features used to test for the differentially expressed genes. Choose from \code{"prev"} (zero/non-zero proportions), \code{"nzm"} (mean), \code{"nzsd"} (non-zero standard deviation), \code{"sd"} (overall standard deviation). Default: \code{"prev"}, \code{"nzm"} and \code{"nzsd"}.}
  \item{verbose}{
Logical. Should the function print the processes? Default: \code{TRUE}.}
  \item{sequencing.data}{
Logical. Is the data.mat a sequencing data matrix (e.g. scRNA sequencing data/count data)? Default: \code{TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{call}{How was the function called?}
\item{R2}{Description of R2}
\item{F0}{Description of F0}
\item{RSS}{Description of RSS}
\item{df.model}{Description of df.model}
\item{df.residual}{Description of df.residual}
\item{coef.list}{Description of coef.list}
\item{p.raw}{Raw, unadjusted P-values.}
\item{p.adj.fdr}{P-values which have been adjusted for false discovery rate.}
\item{p.adj.fwer}{P-values which have been adjusted for family-wise error rate.}
}
\references{
Zhang, L., Yang, L., Ren, Y., Zhang, S., Guan, W., & Chen, J. (Submitted to Bioinformatics): DiSC: a Statistical Tool for Fast Differential Expression Analysis of Individual-level Single-cell RNA-seq Data.
}
\author{
Jun Chen <\email{chen.jun2@mayo.edu}> and Lujun Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
set.seed(seed = 1234556)
data(sim_data)

count_matrix <- sim_data$count_matrix
meta_cell <- sim_data$meta_cell
gene_index <- sim_data$gene_index
meta_ind <- sim_data$meta_ind

obj1 <- DiSC(data.mat = count_matrix, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 999,
             features = c('prev', 'nzm', 'nzsd'), verbose = TRUE,
             sequencing.data = TRUE)
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

obj2 <- DiSC(ct.mat = rarefy_mat, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 99,
             features = c('prev', 'nzm', 'nzsd'), verbose = FALSE)
obj1$p.adj.fdr[123]
obj2$p.adj.fdr[123]
sum(obj1$p.adj.fdr[gene_index$mean_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$var_index] <= 0.10)
sum(obj1$p.adj.fdr[gene_index$mean_var_index] <= 0.10)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
