\name{sim_data}
\alias{sim_data}
\docType{data}
\title{
Simulated Single-cell Sequencing Dataset
}
\description{
This dataset was generated in the "Generate Simulation Datasets" step in the Parametric_simulation.rmd(\url{https://github.com/Lujun995/DiSC_Replication_Code})
}
\usage{data("sim_data")}
\format{
  It contains 12 cases and 12 controls, each with 375 cell replicates. The read depths of each cell replicate are well-balanced. A covariate called RIN (RNA Integrity Number) at the individual level is included in the dataset.

  The dataset comprises a total of 8,000 genes. The signal density was 15\%, with differences in mean, variance, and mean+variance (each at 5\%). The ground truth of differential/equally expression genes are indicated by \code{gene_index}, including \code{mean_index} (genes with a difference in mean), \code{var_index} (genes with a difference in variance), \code{mean_var_index} (genes with a difference in both mean and variance), \code{EE_index} (otherwise (to estimate type-I error)).

  \describe{
  A list of elements:
    \item{\code{count_matrix}}{A numeric count matrix.}
    \item{\code{meta_cell}}{A data.frame of the metadata at the cell level.}
    \item{\code{meta_ind}}{A data.frame of the metadata at the individual level.}
    \item{\code{gene_index}}{A list of 4 numeric vectors representing the ground truth of the IDs of the differentially or equally expressed genes.}
    \item{\code{random_index_gene}}{Not used.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Simulated in the "Generate Simulation Datasets" step in the Parametric_simulation.rmd(\url{https://github.com/Lujun995/DiSC_Replication_Code})
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(sim_data)
str(sim_data)
}
\keyword{datasets}
