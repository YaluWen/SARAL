#' Example dataset for SARAL package
#'
#' This is an example of how to use SARAL to build risk prediction models.
#' Data is obtained from Alzheimer's Disease Neuroimaging Initiatives (ADNI).
#' This example only contains partial data from the ADNI study
#' and thus the result from this analysis should not be interpreted
#' as the results from the original study.
#'
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{Geno}{Genotypes of all individuals (including both the training and testing data).}
#'   \item{Genelist}{A list of 96 genes. Each gene consists a number of single nucleotide variant.}
#'   \item{TrainID}{A bool vector indicating whether the observation is in the training sample.}
#'   \item{Y}{Phenotypes of each individual.}
#' }
#' @source \url{http://adni.loni.usc.edu}
"example_SARAL"
