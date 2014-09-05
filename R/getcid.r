#' Cluster IDs of probes/genes from fitted \code{\linkS4class{Rnits}}
#' 
#' Retrieve cluster IDs of probes/genes from fitted \code{\linkS4class{Rnits}} object after \code{fit} has been run.
#' 
#' If \code{cluster = False} during fitting, a vector of \code{1s} will be returned.
#' 
#' @export
#' @docType methods
#' @rdname getCID-methods
#' @return A vector of cluster IDs corresponding to gene/probe names
#' @param object \code{\linkS4class{Rnits}}
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#' 
#' # Get cluster IDs from fitted model
#' cid <- getCID(rnitsobj)
#'}
setGeneric("getCID", function(object) {
    standardGeneric("getCID")
})
#' @rdname getCID-methods
#' @aliases getCID,character,ANY-method
setMethod("getCID", signature = "Rnits", function(object) {
    x <- object@fitData$fit$clusterID
    x <- setNames(x, rownames(object@fitData$fit))
    return(x)
}) 
