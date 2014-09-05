#' Calculate the optimal B-spline model using generalized cross-validation
#' 
#' Calculate the optimal B-spline model using generalized cross-validation
#' 
#' The optimal B-spline model is chosen as the largest model that minimizes the cross validation error of the top N eigenvectors of each time series data. 
#' 
#' @export
#' @docType methods
#' @rdname calculateGCV-methods
#' @param object \code{\linkS4class{Rnits}} object
#' @param topcomp The number of top eigenvectors to be used for computation
#' 
#' @return A list object with fields 'degree', 'df' for each time series data set.
setGeneric("calculateGCV", function(object, topcomp = 5) {
    standardGeneric("calculateGCV")
})
#' @rdname calculateGCV-methods
#' @aliases calculateGCV,character,ANY-method
#' @importFrom splines bs
#' @importFrom boot cv.glm
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \donttest{
#' opt_model <- calculateGCV(rnitsobj)
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, model = opt.model)
#' }

setMethod("calculateGCV", signature = "Rnits", function(object, topcomp = 5) {
    callData <- object@callData
    Nc <- callData$Nc
    Nsets <- callData$Nsets
    samples <- callData$sample_names
    breaks <- callData$breaks
    timevec <- callData$time_vec
    control <- callData$control_samp
    # Construct matrix of model designs
    maxdf <- min(Nc - 2)
    dfvec = degvec = c()
    for (df in seq(3, maxdf)) {
        dfvec = c(dfvec, df:maxdf)
        degvec = c(degvec, rep(df - 1, length(df:maxdf)))
    }
    idlist = data.frame(deg = degvec, df = dfvec)
    ## SVD
    lr <- exprs(object)
    ## Impute data. Not necessary for time-aligned.
    sinkfile = "/dev/null"
    if (.Platform$OS.type == "windows") {
        sinkfile = "NUL"
    }
    sink(sinkfile)
    if (any(is.na(lr))) 
        lr <- suppressWarnings(impute.knn(lr)$data)
    sink()
    breaks <- c(breaks, ncol(lr) + 1)
    optim_model <- list()
    for (n in 1:(Nsets + 1)) {
        sample_data <- lr[, breaks[n]:(breaks[n + 1] - 1)]
        sample_svd <- svd(t(scale(t(sample_data), scale = FALSE)))
        eigenvalue <- sample_svd$d^2/sum(sample_svd$d^2)
        sig_evs <- which(cumsum(eigenvalue) >= 0.8)[1]
        pred_mat <- c()
        for (vec in 1:topcomp) {
            pred_err <- c()
            for (id in 1:nrow(idlist)) {
                X = bs(timevec, degree = idlist[id, 1], df = idlist[id, 2], knots = placeKnots(sample_data, 
                  idlist[id, 1], idlist[id, 2], timevec), intercept = TRUE)
                temp_dat <- data.frame(X, evec = sample_svd$v[, vec])
                fitglm <- glm(evec ~ ., family = gaussian, temp_dat)
                cvglm <- cv.glm(temp_dat, fitglm)
                pred_err <- c(pred_err, cvglm$delta[2])
            }
            pred_mat <- rbind(pred_mat, pred_err)
        }
        comp_min_errors <- apply(pred_mat, 2, function(x) min(x))
        min_err <- which.min(comp_min_errors)
        optim_model[[paste0("Sample_", n)]] <- data.frame(degree = idlist[min_err, 
            1], df = idlist[min_err, 2], pred_error = min(comp_min_errors), eig1 = eigenvalue[1], 
            sig_evs = sig_evs)
    }
    return(optim_model)
}) 
