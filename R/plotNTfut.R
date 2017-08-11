#' Plotting N-T Plane for Decision on Performing an Interim Analysis
#'
#' The package plots N-T plane for decision for conducting an interim analysis in a randomized control trial.
#' The functions for interim analysis expecting early stopping for superiority and futility are prepared respectively.
#'
#'
#' @usage
#' plotNT.fut(p0, M, q, alpha1, cp1,
#'             xlab = "N: Number of patients at interim analysis",
#'             ylab = "T: Number of responders at interim analysis",
#'             col = "blue",
#'             main = "N-T plot",
#'             lty = 2,...)
#'
#' @param p0              Expected response rate for the control arm: scalar or vector. the value between 0 to 1.
#' @param M               Total number of patients: expected number of patients until last time.
#' @param q               Allocation ratio of the treatment arm: the value between 0 to 1.
#' @param alpha1          Critical alpha at an interim analysis.
#' @param cp1             Critical conditional power at an interim analysis.
#' @param xlab            Label name for x-axis in N-T plot.
#' @param ylab            Label name for y-axis in N-T plot.
#' @param col             Line color. Default is "blue". For multiple p0, set the same length of color with p0.
#' @param main            Main title in N-T plot.
#' @param lty             Line type. The default is 2 for early stopping for futility.
#' @param ...             Other graphics parameters
#'
#' @return A matrix or list with variable names \code{N}, \code{T}, \code{Z_score} and \code{CP}.
#' @return Draw N-T plot for early stopping for futility
#' @return \item{x axis:}{N (Total number of patients at interim analysis)}
#' @return \item{y axis:}{T (Total number of responders at interim analysis)}
#' @details For more details, please refer to the vignette: \code{browseVignettes(package = "IAbin")}
#'
#' @importFrom graphics legend par plot
#' @importFrom stats pnorm qnorm
#' @references
#' Decision on Performing Interim Analysis for Comparative Clinical Trials
#'
#'
#' @examples
#' #--- Settings ---#
#' #--- With an expected parameter for control therapy ---#
#' p0 = 0.5
#' M = 135
#' q = 2/3
#' alpha1 = 0.01
#' cp1 = 0.2
#'
#' #--- N-T plot for early stopping for superiority and futility ---#
#' NT_f = plotNT.fut(p0, M, q, alpha1, cp1)
#' print(NT_f)
#'
#' #--- Settings ---#
#' #--- With several expected parameters for control therapy ---#
#' p0 = c(0.2, 0.4, 0.6)
#' M = 135
#' q = 2/3
#' alpha1 = 0.01
#' col = c(1, 2, 3)
#' cp1 = 0.2
#'
#' #--- N-T plot for early stopping for superiority and futility ---#
#' NT_f3 = plotNT.fut(p0, M, q, alpha1, cp1, col=col)
#' print(NT_f3)
#'
##===============================================================##
##  Function for plotting the critical N-T for the expected p0   ##
##              For early stopping for Futility                  ##
##===============================================================##
#' @export
plotNT.fut = function(p0, M, q, alpha1, cp1,
                      xlab = "N: Number of patients at interim analysis",
                      ylab = "T: Number of responders at interim analysis",
                      col = "blue",
                      main = "N-T plot",
                      lty = 2,...){

  if(is.null(p0)){
    stop("Please set any expected p0 for conrol arm.")
  }
  if(length(p0) != length(col)){
    stop("Please set the same length for the color vectors with the length of p0.")
  }


  if(!is.null(p0)){

    if(length(p0) == 1){
      obj = NT_plane_fut(p0, M, q, alpha1, cp1)
      plot(obj[,1], obj[,2], lty = lty, type = "l",
           xlab = xlab, ylab = ylab,
           xlim = c(0, M), ylim = c(0, M), col = col,
           main = main, cex.main = 1.5, lwd = 2)

      legend("topleft", legend = c(sprintf("%.2f", p0[1])), lty = lty,
             title = expression(paste("Size of ", p[0])), col = col, lwd = 3, bty = "n")

    }else{

      obj = list(NA)
      step = length(p0)

      obj[[1]] = NT_plane_fut(p0[1], M, q, alpha1, cp1)
      plot(obj[[1]][,1], obj[[1]][,2], lty = lty, type = "l",
           xlab = xlab, ylab = ylab,
           xlim = c(0, M), ylim = c(0, M), col = col[1],
           main = main, cex.main = 1.5, lwd = 3)

      for(i in 2:step){
        obj[[i]] = NT_plane_fut(p0[i], M, q, alpha1, cp1)
        par(new = T)
        plot(obj[[i]][,1], obj[[i]][,2], lty = lty, type = "l",
             xlab = "", ylab = "",
             xlim = c(0, M), ylim = c(0, M), col = col[i],
             main = "", lwd = 3)
      }

      leg.val = NULL
      for(i in 1:step){
        leg.val[i] = sprintf("%.2f", p0[i])
      }
      legend("topleft", legend = leg.val, lty = lty,
             title = expression(paste("Size of ", p[0])), col = col, lwd = rep(step, 3), bty = "n")
    }
    options(warn = -1)
    return(invisible(obj))
  }
}




##=================================================================##
##  Function for calculating the critical N-T for the expected p0  ##
##              For early stopping for Futility                    ##
##               Binomial ver. for response rate        (hidden)   ##
##=================================================================##
NT_plane_fut = function(p0, M, q, alpha1, cp1){

  ## Initial check.
  if(q <= 0 || q >= 1){
    stop("Please specify q, the allocation ratio of the new treatment, from 0 to 1.")
  }

  M0 = M*(1-q)
  z_alpha = qnorm(1 - alpha1/2, 0, 1)
  z_gamma = qnorm(1 - cp1, 0, 1)

  ## Calculate the expected z score for each T.
  reslist = list(NA)

  for(N0 in 1:M0){
    res = matrix(NA, nrow= round(N0 / (1 - q)), ncol = 4)
    colnames(res) = c("T", "deltahat", "z_stat", "z_fut")
    for(T in 1:round(N0 / (1 - q))){
      N1 = round(q / (1 - q) * N0)
      N01 = round(N0 / (1 - q))
      p1hat = (T - N0 * p0) / (N1)
      deltahat = p1hat - p0
      Rhat = q * p1hat + (1 - q) * p0
      var_p1hat = (N01 * Rhat * (1 - Rhat)) / ((N01 * q)^2)
      z_stat = deltahat / sqrt(var_p1hat)
      z_fut = z_alpha * sqrt(M / N01) - z_gamma * sqrt((M - N01) / N01) - (deltahat * (M - N01)) / sqrt(var_p1hat * N01)
      res[T, ] = c(T, deltahat, z_stat, z_fut)
    }
    reslist[[N0]] = res
  }

  ## search for "smallest T" such that z score <= z_fut.
  res_plane = matrix(NA, nrow = M0, ncol = 4)
  colnames(res_plane) = c("N", "T", "Z_score", "CP")

  for(i in 1:M0){
    mat = as.data.frame(reslist[[i]])
    vec = mat[mat$z_stat <= mat$z_fut, ]
    Tmax = max(vec$T, na.rm = T)
    N01 = round(i/(1 - q))
    z_score = max(vec$z_stat, na.rm = T)
    deltahat = max(vec$deltahat, na.rm = T)
    phat = Tmax / N01
    varhat1 = phat * (1 - phat) * (1 / i + 1 / (N01 - i))
    c = (z_alpha * sqrt(varhat1) * sqrt(M) - deltahat * sqrt(N01)) / (M - N01)
    zg = (c - deltahat) / (sqrt(varhat1) / sqrt(M - N01))
    CP = 1 - pnorm(zg, 0, 1)
    res_plane[i, ] = c(N = N01, Tmax, z_score, CP)
  }
  options(warn = -1)
  return(res_plane)
}
