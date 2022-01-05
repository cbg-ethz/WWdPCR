#' Swiss clinical data
#'
#' A dataset containing Swiss clinical sequencing data used in the publication:
#' Caduff, Lea, David Dreifuss, Tobias Schindler, Alexander J. Devaux,
#' Pravin Ganesanandamoorthy, Anina Kull, Elyse Stachler et al.
#' "Inferring Transmission Fitness Advantage of SARS-CoV-2 Variants of Concern
#' in Wastewater Using Digital PCR."
#'  medRxiv (2021).
#'
#' @format A data frame with 4037 obs. of  6 variables.
#' \describe{
#'   \item{division}{factor with 27 levels, the Swiss canton}
#'   \item{date}{factor with 424 levels, representing the date of sampling}
#'   \item{HV6970del}{the count of sequences with the S-gene HV69-70 deletion}
#'   \item{ORF1adel}{the count of sequences with the ORF1a deletion}
#'   \item{B117}{the count of sequences determined to be B.1.1.7 (Alpha)}
#'   \item{Tot}{the total count of sequences}
#' }
#' @source \url{https://cov-spectrum.org}
"clinical_data"


#' Zurich HV69-70 dPCR data
#'
#' A dataset containing Zurich HV69-70 dPCR data used in the publication:
#' Caduff, Lea, David Dreifuss, Tobias Schindler, Alexander J. Devaux,
#' Pravin Ganesanandamoorthy, Anina Kull, Elyse Stachler et al.
#' "Inferring Transmission Fitness Advantage of SARS-CoV-2 Variants of Concern
#' in Wastewater Using Digital PCR."
#'  medRxiv (2021).
#'
#' @format A data frame with 64 obs. of  6 variables.
#' \describe{
#'   \item{date}{Date, representing the date of sampling}
#'   \item{dayssince}{number of days since 2020-12-07}
#'   \item{Dilutionfactor}{factor used in dilution,
#'   useful for calculating viral RNA concentration}
#'   \item{x0}{double negative droplets count}
#'   \item{x1}{single positive droplets count}
#'   \item{x2}{double positive droplets count}
#' }
#' @source \doi{10.1101/2021.08.22.21262024}
"dPCRdat_HV6970"
