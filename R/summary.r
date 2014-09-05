#' Summary of fit
#' 
#' Summarize top genes or probes from Rnits \code{fit} method
#' 
#' @param object \code{\linkS4class{Rnits}} object on which \code{fit} has been applied
#' @param top Display results for top N genes/probes. Default \code{50}
#' @param fdr Display results for genes/probes less than FDR (\%) cutoff (if provided). Overrides \code{top} argument
#' @param plot If \code{TRUE}, plot histogram of p-values
#' @param sort.by Sort top results by either \code{p-value} or \code{FDR}
#' @export
#' @rdname summary-methods
#' @docType methods
#' @return A table of top genes/profiles listing the ratio statistics, p-values, q-values and  cluster information.
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'
#' # Get summary of top genes
#' summary(rnitsobj, FDR = 5)
#'
#'}
setMethod("summary", signature = "Rnits", function(object, top = 48, fdr = NULL, 
    plot = FALSE, sort.by = "p-value") {
    # Check object
    if (class(object) != "Rnits") 
        stop("Object not of Rnits class")
    # Check top
    if (!is.numeric(top)) 
        stop("top must be an integer")
    if (top < 1 | top > 1e+05) 
        top = 48
    # Check fdr
    if (!is.null(fdr)) {
        if (fdr > 100 | fdr < 0) 
            fdr = NULL
    }
    ## Check sort.by
    if (sort.by != "p-value" & sort.by != "FDR") {
        stop("sort.by argument can only take value p-value or FDR")
    }
    ## Initialize
    top <- 1:top
    if (!object@callData$fit.model) 
        stop("Run fit() on data first")
    gene.level <- object@callData$gene.level
    if (gene.level) {
        name <- as.vector(unique(fData(object)$GeneName))
    } else {
        name <- featureNames(object)
    }
    ## Get Stats
    pval <- getPval(object)
    ratio <- getStat(object)
    qval <- qvalue(pval)$q
    clusModel <- object@fitData$clusters
    clusID <- object@fitData$fit$clusterID
    ## Construct table
    out <- data.frame(Rank = 1:length(name), Name = name, Statistic = round(ratio, 
        3), p.value = round(pval, 6), `FDR percent` = round(qval * 100, 3), ClusterID = clusID, 
        Model.degree = clusModel[clusID, 1], Model.df = clusModel[clusID, 2])
    colnames(out)[2] <- "Name"
    ## Sort
    if (sort.by == "FDR") {
        out <- out[with(out, order(FDR.percent, -Statistic)), ]
    } else out <- out[with(out, order(p.value, -Statistic)), ]
    ## Plot histogram
    if (plot) 
        hist(pval, main = "Histogram of p-values", xlab = "p-values", col = "gray")
    if (!is.null(fdr)) 
        top <- which(out$FDR.percent <= fdr)
    out.trunc <- out[top, ]
    out.trunc$Rank <- 1:length(top)
    return(out.trunc)
}) 
