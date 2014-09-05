#' Data of top genes/probes
#' 
#' Extract expression data for top genes/probes
#' 
#' @param object \code{\linkS4class{Rnits}} object on which \code{fit} has been applied
#' @param id Names of probes or genes
#' @param top Display results for top N genes/probes. Default \code{50}
#' @param fdr Display results for genes/probes less than FDR cutoff (if provided). Overrides \code{top} argument
#' @param sort.by Sort top results by either \code{p-value} or \code{FDR}
#' 
#' @return A table of expression values of top genes/profiles 
#' @export
#' @docType methods
#' @rdname topData-methods
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'
#' #Get data for top genes
#' td <- topData(rnitsobj, FDR = 5)
#' 
#'}
#'
#'
setGeneric("topData", function(object, id = NULL, fdr = NULL, top = 16, sort.by = "p-value") {
    standardGeneric("topData")
})
#' @rdname topData-methods
#' @aliases topData,character,ANY-method
setMethod("topData", signature = "Rnits", function(object, id = NULL, fdr = NULL, 
    top = 16, sort.by = "p-value") {
    if (!object@callData$fit.model) 
        stop("Run fit() on data first")
    lr <- getLR(object)
    genes <- rownames(lr)
    timevec = sort(unique(pData(object)$Time))
    colnames(lr) <- paste0(pData(object)$Sample, "_", timevec)
    if (is.null(fdr)) {
        stable <- summary(object, top = top, sort.by = sort.by)
    } else {
        stable <- summary(object, fdr = fdr, sort.by = sort.by)
        cat(nrow(stable), "\n")
    }
    topNames <- as.vector(stable$Name)
    hits <- match(toupper(topNames), toupper(genes))
    return(lr[hits, ])
}) 
