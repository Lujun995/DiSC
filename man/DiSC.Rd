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
DiSC(ct.mat, cell.ind, metadata, outcome, covariates = NULL,
     stats.combine.func = max, cell.id = "cell_id",
     individual.id = "individual", perm.no = 999,
     features = c('prev', 'nzm', 'nzsd'), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct.mat}{
Sequencing count matrix. Genes are in rows, column names should be cell ids.
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
  \item{stats.combine.func}{
A function defines the method to combine F-statistics using different features, default \code{max}.}
  \item{cell.id}{
A character string which is the variable name of cell ids in \code{cell.ind}.}
  \item{individual.id}{
A character string which is the variable name of the individual id variables in \code{cell.ind} and \code{metadata}.}
  \item{perm.no}{
An integer, number of permutations used. Default: \code{999}. It can be reduced to \code{99} if adjusted P-values (false discovery rate) are the only interest.}
  \item{features}{
Vector of a character string. Features used to test for the differentially expressed genes. Choose from \code{"prev"} (zero/non-zero proportions), \code{"nzm"} (non-zero mean), \code{"nzsd"} (non-zero standard deviation), \code{"sd"} (overall standard deviation). Default: \code{"prev"}, \code{"nzm"} and \code{"nzsd"}.}
  \item{verbose}{
Logical. Should the function print the processes? Default: \code{TRUE}.}
  \item{...}{
Other arguments passed to function \code{basic.func}.}
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
%% ~put references to the literature/web site here ~
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
\code{\link[GUniFrac]{Rarefy}}
}
\examples{
\dontrun{
set.seed(seed = 1234556)
data(sim_data)

count_matrix <- sim_data$count_matrix
meta_cell <- sim_data$meta_cell
gene_index <- sim_data$gene_index
meta_ind <- sim_data$meta_ind

rarefy_mat <- t(GUniFrac::Rarefy(t(count_matrix))$otu.tab.rff)

obj1 <- DiSC(ct.mat = rarefy_mat, cell.ind = meta_cell,
             metadata = meta_ind, outcome = "phenotype",
             covariates = "RIN", cell.id = "cell_id",
             individual.id = "individual", perm.no = 999,
             features = c('prev', 'nzm', 'nzsd'))
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
