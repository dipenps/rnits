#' Get log-ratios
#' 
#' Extract normalized log-ratios from \code{\linkS4class{Rnits}} object
#' @aliases Rnits.getLR 
#' @param object \code{\linkS4class{Rnits}} object
#' @param impute If \code{TRUE}, perform K-NN imputation to fill missing values
#' @return A matrix of normalized log-ratios.
#' @export
#' @docType methods
#' @rdname getLR-methods
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \dontrun{
#' # Fit model using gene-level summarization
#' rnitsobj <- fit(rnitsobj, gene.level = TRUE, clusterallsamples = FALSE)
#'
#' # Get logratios
#' lr <- getLR(rnitsobj)
#'}
setGeneric("getLR", function(object, impute = FALSE) {
    standardGeneric("getLR")
})
#' @rdname getLR-methods
#' @aliases getLR,character,ANY-method
setMethod("getLR", signature = "Rnits", function(object, impute = FALSE) {
    gene.level = object@callData$gene.level
    if (!is.null(gene.level) && gene.level) {
        lr = exprs(object@geneData)
    } else {
        lr = exprs(object)
    }
    
    if(impute) lr <- impute.knn(as.matrix(lr))$data
    return(lr)
}) 
