#' Fit model on time series data
#' 
#' Fit a model comparing time series data set \code{\linkS4class{Rnits}} objects
#' 
#' The function compares multiple time-series expression data sets by i) (optional) 
#' summarizing probes into gene-level information ii) (optional) identifying a set 
#' of co-expressed genes by clustering iii) For each cluster (or for all genes
#' /probes), fit a series of B-splines with varying curvature and degrees of 
#' freedom. Under the null hypothesis \code{H_0}, a single model is fit for 
#' all data sets, while under \code{H_1}, each data set is fit separately. 
#' P-values from the hypothesis test are then plotted and the least complex 
#' spline parameters that result in uniformly distributed null p-values are 
#' automatically chosen.  
#' @param object \code{\linkS4class{Rnits}} object
#' @param cluster if \code{TRUE}, perform clustering to identify groups of 
#' genes/probes with similar expression profiles.
#' @param B Default \code{100}. Number of bootstrap iterations for p-value 
#' calculation
#' @param verbatim If \code{FALSE}, print out details of fitting models.
#' @param nclus Default \code{NULL}. Number of clusters to use for k-means 
#' clustering.
#' @param modelhistplot If \code{TRUE}, p-value histograms of multiple 
#' models are plotted. 
#' @param seed Random seed for bootstrap iterations
#' @param gene.level If \code{TRUE}, collapse probes to gene level information.
#' @param clusterallsamples If \code{TRUE}, Use all time series for clustering. 
#' By default, only the sample labeled 'control' is used or the lexically 
#' first sample is used.
#' @param model A data frame with fields 'degree' and 'df' indicating a specific 
#' B-spline model to be used. If provided, model selection is not run.
#' 
#' @export
#' @docType methods
#' @rdname fit-methods
#' @return An object of S4 class \code{\linkS4class{Rnits}} with fitted results 
#' data containing cluster information, ratio statistics and p-values.
#' @importFrom qvalue qvalue

#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'}
#'
setGeneric("fit", function(object, cluster = TRUE, B = 100, verbatim = FALSE, nclus = NULL, 
                           modelhistplot = FALSE, seed = 123, gene.level = TRUE, clusterallsamples = FALSE, 
                           model = NULL) {
  standardGeneric("fit")
})

#' @rdname fit-methods
#' @aliases fit,character,ANY-method
setMethod("fit", "Rnits", function(object, cluster = TRUE, B = 100, verbatim = FALSE, 
                                   nclus = NULL, modelhistplot = FALSE, seed = 123, gene.level = TRUE, clusterallsamples = FALSE, 
                                   model = NULL) {
  
  ## Check object
  if (class(object) != "Rnits") 
    stop("Object not of Rnits class")
  if (!is.null(object@callData$fit.model) && object@callData$fit.model) 
    stop("The object has been model-fit before. Use new instance for re-analysis.")
  
  ## Check cluster
  if (!cluster %in% c(TRUE, FALSE)) 
    cluster = TRUE
  
  ## Check B
  if (B < 1 | B > 1000) 
    B = 100
  
  ## Check modelhistplot
  if (!modelhistplot %in% c(TRUE, FALSE)) 
    modelhistplot = FALSE
  
  ## Check gene.level
  if (!is.null(object@callData$gene.level)) 
    gene.level <- object@callData$gene.level
  if (!gene.level %in% c(TRUE, FALSE)) 
    gene.level = FALSE
  
  ## Check if gene level or probe level summarization
  if (gene.level) {
    cat("Collapsing probes into gene level log-ratios\n")
    object <- summarizeProbes(object)
    lr <- exprs(object)
    object@callData$gene.level <- TRUE
  } else {
    object@callData$gene.level <- FALSE
  }
  
  ## If time alignment is not done previously, extract log-ratios from object
  if (!gene.level) {
    lr <- exprs(object)
  } else {
    lr <- exprs(object@geneData)
  }
  ## Impute data. Not necessary for time-aligned.
  set.seed(seed)
  sinkfile = "/dev/null"
  if (.Platform$OS.type == "windows") {
    sinkfile = "NUL"
  }
  sink(sinkfile)
  if (any(is.na(lr))) 
    lr <- suppressWarnings(impute.knn(lr)$data)
  sink()
  
  if (!gene.level) {
    exprs(object) <- lr
  } else {
    exprs(object@geneData) <- lr
  }
  
  
  Ng = nrow(lr)
  
  callData <- object@callData
  Nc <- callData$Nc
  Nrep <- callData$Nrep
  Nsets <- callData$Nsets
  samples <- callData$sample_names
  breaks <- callData$breaks
  timevec <- rep(callData$time_vec, Nrep)
  control <- callData$control_samp
  
  ## Construct matrix of model designs
  maxdf <- min(Nc/Nrep - 2, 8)
  dfvec = degvec = c()
  for (df in seq(3, maxdf)) {
    dfvec = c(dfvec, df:maxdf)
    degvec = c(degvec, rep(df - 1, length(df:maxdf)))
  }
  idlist = data.frame(deg = degvec, df = dfvec)
  if (max(idlist[, 1]) > 5) 
    idlist <- idlist[-which(idlist[, 1] > 5), ]
  
  if (!is.null(model)) {
    cat("Using specified model parameters\n")
    idlist = data.frame(deg = model$degree, df = model$df)
  }
  cat("############################", "\n")
  cat("Analyzing Models", "\n")
  pcomb = ratiocomb = rep(0, Ng)
  
  
  ## Compute Clusters
  clustersamples = control
  if (clusterallsamples == TRUE) 
    clustersamples = 1:ncol(lr)
  
  if (cluster == TRUE) {
    if (is.null(nclus)) {
      csvd = svd(t(scale(t(lr[, clustersamples]), scale = FALSE)))
      d = (csvd$d^2)/sum(csvd$d^2)
      nclus = min(2 * which(cumsum(d) > 0.7)[1], 16)
      
      ## Gap statistic (Depends on package 'cluster') clusgap =
      # clusGap(lr[sample(nrow(lr), 500), clustersamples], FUN = kmeans, K.max = 12, B
      # = 100) nclus = maxSE(clusgap$Tab[,3], clusgap$Tab[,4])
      cat("Computed number of clusters is", nclus, "\n")
    } else {
      cat("Using user-specified number of clusters of ", nclus, "\n")
    }
    k.clus = kmeans(lr[, clustersamples], centers = nclus, nstart = 10)
    while (sum(table(k.clus$cluster) < 500) & nclus > 1) {
      cat("Repeating clustering to eliminate small clusters\n")
      nclus <- max(1, nclus - 1)
      k.clus = kmeans(lr, centers = nclus, nstart = 1)
    }
    k.clus = kmeans(lr, centers = nclus, nstart = 10)
  }
  
  if (cluster == FALSE) {
    k.clus <- list()
    k.clus$cluster = rep(1, Ng)
    nclus = 1
  }
  cat(paste("Number of clusters is:", nclus), "\n")
  cat("############################", "\n")
  ##################### Begin analysis
  comb = ratiocomb = rep(0, Ng)
  clus.models = matrix(0, nclus, 2)
  for (nc in 1:nclus) {
    id.nc = which(k.clus$cluster == nc)
    cat("Cluster", nc, "of size", length(id.nc), "\n")
    inpdata = lr[id.nc, ]
    pcomb.clus = rcomb.clus = p.ks.clus = pi0.clus = c()
    for (id in 1:nrow(idlist)) {
      if (verbatim) 
        cat("Model", id, "|", "Degree =", idlist[id, 1], ", DF = ", idlist[id, 2], "\n") 
      
      
      ## Set basis function
      X = bs(timevec, degree = idlist[id, 1], df = idlist[id, 2], 
             knots = placeKnots(inpdata[, control], idlist[id, 1], idlist[id, 2], timevec), intercept = TRUE) 
      
      
      ## Fit model
      s0 = 0  # Initially compute s0
      modelFit = tsFit(inpdata, X, breaks, s0, verbatim = FALSE, B = B)
      ratio = modelFit$ratio
      ratiostar = modelFit$ratiostar
      pcomb.clus = cbind(pcomb.clus, getp(ratio, ratiostar))
      rcomb.clus = cbind(rcomb.clus, ratio)
      
      # Bioc version 3.9. pi0 estimation is deprecated (March 2019)
      #pi0.clus <- c(pi0.clus, qvalue(getp(ratio, ratiostar))$pi0)
      pi0.clus <- c(pi0.clus, 1)
      if (verbatim) 
        cat(" || computed s0 is ", round(modelFit$s0, digits = 3), "\n")
      if (verbatim) 
        cat(" ||| estimated pi0 is ", pi0.clus[length(pi0.clus)], "\n")
      
      ##################### Calculate KS-test for uniform distribution of p-values ix =
      ##################### which(pcomb.clus[,id] >= quantile(pcomb.clus[,id], .75)) if(sum(pcomb.clus[,id]
      ##################### > 0.5) > 100) { tH = hist(pcomb.clus[which(pcomb.clus[,id]>= 0.5),id], plot =
      ##################### FALSE, breaks = 20) p.ks.clus[id] =
      ##################### summary.lm(lm(tH$density~tH$mids))$coef[2,1] } else if(sum(pcomb.clus[,id] >
      ##################### 0.25) > 100){ tH = hist(pcomb.clus[which(pcomb.clus[,id]>= 0.20),id], plot =
      ##################### FALSE, breaks = 20) p.ks.clus[id] =
      ##################### summary.lm(lm(tH$density~tH$mids))$coef[2,1] } else { if(verbatim) cat('Low
      ##################### p-values. Model Selection may be unstable.\n') p.ks.clus[id] = 10 } }
      ##################### if(all(p.ks.clus == 10)) { cat('Low p-values for all candidate models. Smallest
      ##################### model selected.\n') mod.pick = 1 } else { mod.pick = which.min(abs(p.ks.clus))
      ##################### 
    } 
    # Select model based on power
    siggene_count <- apply(pcomb.clus, 2, function(x) sum(p.adjust(x, "fdr") >= 0.1))
    
    if (sum(siggene_count == min(siggene_count)) == 1) {
      mod.pick <- which.min(siggene_count)
    } else {
      mod.pick <- which(siggene_count == min(siggene_count))[1]
    }
    
    if (verbatim) 
      cat("Selected model is Model ", mod.pick, " | Degree = ", idlist[mod.pick, 1], " DF = ", idlist[mod.pick, 2], "\n")
    if (modelhistplot) {
      par(mfrow = c(3, 2))
      for (n in 1:ncol(pcomb.clus)) {
        hist(pcomb.clus[, n], 50, xlab = "pvalues", main = paste0("Model ", n))
      }
    }
    
    # mod.pick = which.min(abs(p.ks.clus))
    if (verbatim) 
      cat("Selected model for Cluster", nc, "is", paste(idlist[mod.pick, 1], "_", idlist[mod.pick, 2], sep = ""), "\n")
    if (verbatim) 
      cat("############################", "\n")
    
    clus.models[nc, ] = c(idlist[mod.pick, 1], idlist[mod.pick, 2])
    pcomb[id.nc] = pcomb.clus[, mod.pick]
    ratiocomb[id.nc] = rcomb.clus[, mod.pick]
  }
  
  #cat("Estimated proportion of null genes is", round(100 * qvalue(pcomb)$pi0, 2), "%\n")
  cat("Estimated proportion of null genes is", round(100 * 1, 2), "%\n")

  object@fitData$fit <- data.frame(Ratio.statistic = ratiocomb, p.value = pcomb, 
                                   clusterID = k.clus$cluster, row.names = rownames(lr))
  object@fitData$clusters <- data.frame(clusterModel.degree = clus.models[, 1], clusterModel.df = clus.models[, 2])
  
  object@callData$fit.model <- TRUE
  
  
  return(object)
}
) 
