#' Get p-values
#' 
#' Extract p-values from fitted \code{\linkS4class{Rnits}} object
#' 
#' @param object \code{\linkS4class{Rnits}} object
#' @return An vector of p-values
#' @export
#' @docType methods
#' @rdname getPval-methods
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'
#' #Get pvalues from fitted model
#' pval <- getPval(rnitsobj)
#' 
#'}
setGeneric("getPval", function(object) {
    standardGeneric("getPval")
})
#' @rdname getPval-methods
#' @aliases getPval,character,ANY-method
setMethod("getPval", signature = "Rnits", function(object) {
    x <- object@fitData$fit$p.value
    x <- setNames(x, rownames(object@fitData$fit))
    return(x)
}) 
