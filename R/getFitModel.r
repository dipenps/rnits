#' Extract fit data from \code{\linkS4class{Rnits}} object
#' 
#' Retrieve model fit data from \code{\linkS4class{Rnits}} object after \code{fit} has been run.
#' 
#' Contains Ratio statistic, p-value and cluster ID data
#' 
#' @export
#' @docType methods
#' @rdname getFitModel-methods
#' @return A data frame containing the model fit results for all genes
#' @param object \code{\linkS4class{Rnits}}
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'
#' # P-values, ratio statistics and cluster ID's can be retrieved for all genes together
#' fitdata <- getFitModel(rnitsobj)
#'}
setGeneric("getFitModel", function(object) {
    standardGeneric("getFitModel")
})
#' @rdname getFitModel-methods
#' @aliases getFitModel,character,ANY-method
setMethod("getFitModel", signature = "Rnits", function(object) object@fitData$fit) 
