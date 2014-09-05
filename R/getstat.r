#' Retrieve ratio statistics
#' 
#' Extract ratio statistics from fitted \code{\linkS4class{Rnits}} object
#' 
#' @param object \code{\linkS4class{Rnits}} object
#' @return An vector of ratio statistics
#' @export
#' @rdname getStat-methods
#' @docType methods
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#' 
#' # Get ratio statistics from fitted model
#' stat <- getStat(rnitsobj)
#'}
setGeneric("getStat", function(object) {
    standardGeneric("getStat")
})
#' @rdname getStat-methods
#' @aliases getStat,character,ANY-method
setMethod("getStat", signature = "Rnits", function(object) {
    x <- object@fitData$fit$Ratio.statistic
    x <- setNames(x, rownames(object@fitData$fit))
    return(x)
}) 
