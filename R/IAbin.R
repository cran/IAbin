#'  @name IAbin-Package
#'  @alias IAbin-Package
#'  @docType package
#'  @title Plotting N-T Plane for Decision on Performing an Interim Analysis
#'  @description
#'  In randomized-controlled trials, interim analyses are often planned for possible early trial termination to claim superiority or futility of a new therapy. Blinded data also have information about the potential treatment difference between the groups. We developed a blinded data monitoring tool that enables investigators to predict whether they observe such an unblinded interim analysis results that supports early termination of the trial. Investigators may skip some of the planned interim analyses if an early termination is unlikely.
#'  This tool will provide reference information about N: Sample size at interim analysis, and T: Total number of responders at interim analysis for decision on performing an interim analysis.
#'  @author Kyongsun Pak
#'  Maintainer: Kyongsun Pak <pakk@pharm.kitasato-u.ac.jp>
#'
#'  @details
#'  \tabular{ll}{
#'  Package: \tab IAbin\cr
#'  Type: \tab Package\cr
#'  Version: \tab 1.0\cr
#'  Date: \tab 2017-01-05\cr
#'  License: \tab GPL-2\cr
#'  Please check the vignette for details: \code{browseVignettes(package = "IAbin")}
#'  }
#'  @references Decision on Performing Interim Analysis for Comparative Clinical Trials
#'
#'
#'  @examples
#'  #--- Settings for expected trial design ---#
#'  p0 = c(0.2, 0.4, 0.6)
#'  M = 135
#'  q = 2/3
#'  alpha1 = 0.01
#'  cp1 = 0.2
#'  col = c(1, 2, 3)
#'
#'  #--- N-T plot for early stopping for superiority ---#
#'  plotNT.sup(p0, M, q, alpha1, col=col)
#'
#'  #--- Adding N-T plot early stopping for futility ---#
#'  par(new = T)
#'  plotNT.fut(p0, M, q, alpha1, cp1, col=col)
#'
