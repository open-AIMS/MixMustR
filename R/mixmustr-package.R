#' The 'mixmustr' package.
#'
#' @description This is a flexible Bayesian mixture model package written in
#' the probabilistic programming language Stan (<https://mc-stan.org/>) for R
#' (<https://www.r-project.org/>). It estimates source mixing proportions by
#' incorporating simultaneous likelihood evaluation from two independent data
#' streams collected from the mixture of interest: one obtained from chemical
#' tracers/biomarkers (i.e., a single tracer measurement per observation, e.g.,
#' from stable isotopes and fatty acids), and another yielding source
#' composition (e.g., based on eDNA). `mixmustr` also allows for the estimation
#' of an additional, unsampled source component to partially relax the
#' assumption that the mixing proportions from all samples sources should sum
#' up to 1. `mixmustr` should have wide applicability in ecological studies,
#' particularly given the growing usage and availability of eDNA and
#' metabarcoding methods to understand carbon source-sink dynamics as well as
#' trophic interactions.
#' @name mixmustr-package
#' @aliases mixmustr
"_PACKAGE"
