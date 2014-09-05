#' @title
#' Curve registration of time series curves
#' 
#' @description
#' Align multiple time series to the average seris
#' 
#' @param object \code{rnits} object
#' @param iterMax Maximum iterations to be performed
#' @param seed Random seed
#' @param null.frac Fraction of genes that are considered as null
#' @param anchor Sample to be considerded as base for aligning time series. If not provided, the average is used
#' @param rerun If \code{TRUE}, re-align previously aligned data
#' @param plot If \code{TRUE}, plot results
#' 
#' @export
#' @docType methods
#' @rdname timeAlign-methods
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' 
#' # Do curve-registration on data
#' rnitsobj <- timeAlign(rnitsobj)
#'
#'
#'
setGeneric("timeAlign", function(object, iterMax = 5, seed = 123, null.frac = 0.75, 
    anchor = NULL, rerun = FALSE, plot = FALSE) {
    standardGeneric("timeAlign")
})
#' @rdname timeAlign-methods
#' @aliases timeAlign,character,ANY-method
setMethod("timeAlign", signature = "Rnits", function(object, iterMax = 5, seed = 123, 
    null.frac = 0.75, anchor = NULL, rerun = FALSE, plot = FALSE) {
    ## Check object
    if (class(object) != "Rnits") 
        stop("Object not of Rnits class")
    ## Check iterMax_
    if (iterMax > 10 | iterMax < 1) 
        iterMax = 5
    ## Check null.frac
    if (null.frac >= 1 | null.frac <= 0) 
        null.frac = 0.75
    if (!is.null(object@callData$time.align)) {
        cat("The data is previously time-aligned. Set rerun = TRUE to rerun.\n")
        if (!rerun) 
            return(object)
    }
    ## Extract normalized log-ratio data
    lr0 <- exprs(object)
    ## Impute data, if necessary
    if (sum(dim(lr0)) == 0) {
        set.seed(seed)
        sinkfile = "/dev/null"
        if (.Platform$OS.type == "windows") {
            sinkfile = "NUL"
        }
        sink(sinkfile)
        if (sum(is.na(lr0)) > 0) 
            lr0 <- suppressWarnings(impute.knn(lr0)$data)
        sink()
    }
    samples <- unique(pData(object)$Sample)
    if (!is.null(anchor)) {
        if (!anchor %in% samples | length(anchor) > 1) {
            cat(paste0("Sample ", anchor, " not found. Using sample ", samples[1], 
                " as anchor\n"))
            anchor = samples[1]
        }
    } else {
        anchor = samples[1]
    }
    targets <- setdiff(samples, anchor)
    anchordat <- pData(object)$Sample %in% anchor
    ################# Iteration begins control = grep('control', data.iter@designData$Cy5) treatment
    ################# = grep('treatment', data.iter@designData$Cy5)
    for (i in targets) {
        iter = 1
        theta.sample = 0
        lr <- lr0
        while (iter < (iterMax + 1)) {
            cat(paste0("Sample: ", i), "\n")
            cat("======Iteration Number ", iter, "======\n")
            targetdat <- pData(object)$Sample %in% i
            # Pick random selected null genes null.id <- which(getPval(data.iter) >
            # quantile(getPval(data.iter), null.frac)) ## Type A
            l2.norm <- rowSums((lr[, targetdat] - lr[, anchordat])^2)
            null.id <- which(l2.norm < quantile(l2.norm, (1 - null.frac), na.rm = TRUE))
            if (length(null.id) == 0) 
                stop("NULL probe error")
            lrC = lr[null.id, anchordat]
            # Fit Fast Fourier Transform
            fftS <- t(mvfft(t(lr[null.id, targetdat])))
            ArgS <- Arg(fftS)
            ModS <- Mod(fftS)
            # Compute theta adjustment on putative nulls
            theta = seq(-pi/2, pi/2, length = 101)
            mse = rep(0, 101)
            for (k in 1:101) {
                lrS <- Re(t(mvfft(t(ModS * exp((0 + (0 + (0 + (0 + (0+1i))))) * ArgS) * 
                  exp((0 + (0 + (0 + (0 + (0+1i))))) * theta[k])), inverse = TRUE)/sum(anchordat)))
                diff <- (lrS - lrC)
                diff.sq = rowSums(diff^2)
                mse[k] = mean(diff.sq)
            }
            theta.adj = theta[which.min(mse)]
            theta.sample = theta.sample + theta.adj
            if (plot) {
                plot(theta, mse, pch = 19, cex = 0.9, main = paste0(i, ": Iter ", 
                  iter), ylab = "MSE", xlab = "Theta (radian)")
                abline(v = theta.sample, col = "red")
            }
            if (theta.adj == min(abs(theta))) {
                cat("Computed Adjustment for  round is ", round(abs(theta.adj)/pi * 
                  100, 2), " % of cycle period\n")
                cat(paste0("Exiting time alignment for sample ", i, " after ", iter, 
                  " round(s)\n"))
                lr0[, targetdat] <- lr[, targetdat]
                object@callData$theta.adj[i] <- theta.sample
                break
            }
            cat("Computed Adjustment for round is ", round(abs(theta.adj)/pi * 100, 
                2), " % of cycle period\n")
            ## Adjust full data set
            fftS.full <- t(mvfft(t(lr[, targetdat])))
            ArgS.full <- Arg(fftS.full)
            ModS.full <- Mod(fftS.full)
            lrS.adj <- Re(t(mvfft(t(ModS.full * exp((0 + (0 + (0 + (0 + (0+1i))))) * 
                ArgS.full) * exp((0 + (0 + (0 + (0 + (0+1i))))) * theta.adj)), inverse = TRUE)/sum(anchordat)))
            lr[, targetdat] <- lrS.adj
            ## Find null genes based on L2 norm
            l2.norm = rowSums((lr[, anchordat] - lrS.adj)^2)
            null0 <- which(l2.norm < quantile(l2.norm, (1 - null.frac), na.rm = TRUE))
            ## Check if pre- and post- adjustment nulls match
            overlap.null = length(intersect(null0, null.id))
            cat(overlap.null, length(null0), "\n")
            iter = iter + 1
            if (overlap.null > 0.95 * length(null.id)) {
                cat("###### Converged in ", (iter - 1), "iterations ##########\n")
                cat(paste0("Total adjustment for sample is ", round(abs(theta.sample)/pi * 
                  100, 2), " % of cycle period\n"), "\n")
                iter = iterMax + 2
                object@callData$theta.adj[i] <- theta.sample
                lr0[, targetdat] <- lr[, targetdat]
                break
            }
            ## If max iterations reached, reset data to original data.
            if (iter == iterMax + 1) {
                cat("Time series alignment did not converge. No alignment done\n")
                break
            }
        }
        cat("\n\n\n")
    }
    exprs(object) <- lr0
    object@callData$time.align = TRUE
    return(object)
}) 
