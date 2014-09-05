setClassUnion("RGListOrNull", c("RGList", "NULL"))
setClassUnion("AffyBatchOrNull", c("AffyBatch", "NULL"))
setClassUnion("ExpressionSetOrNull", c("ExpressionSet", "NULL"))
#' rnits class
#' 
#' Class rnits for time series
#' 
#' Some details
#' @export
#' @import methods
#' @importClassesFrom limma RGList
#' @importClassesFrom affy AffyBatch
#' 
setClass("Rnits", contains = "ExpressionSet", representation(rawTwoColor = "RGListOrNull", 
    rawAffy = "AffyBatchOrNull", fitData = "list", callData = "list", geneData = "ExpressionSetOrNull"), 
    prototype = prototype(rawTwoColor = NULL, rawAffy = NULL, geneData = NULL))
setValidity("Rnits", function(object) {
    phenodata <- object@phenoData@data
    if (!toupper("Time") %in% toupper(colnames(phenodata))) 
        stop("Pheno Data must contain a Time Column")
    if (!toupper("Sample") %in% toupper(colnames(phenodata))) 
        stop("Pheno Data must contain a Sample Column")
    if (length(unique(table(phenodata$Time))) != 1) {
        stop("Mismatched time data")
    } else {
        cat(paste0("Data set with ", unique(table(phenodata$Time)), " time series\n"))
    }
    if (!all(rownames(object@featureData@data) == rownames(exprs(object)))) {
        stop("Feature names do not match with assay data")
    }
})
setMethod("initialize", "Rnits", function(.Object, ...) {
    callNextMethod(.Object, ...)
})
#' replace slot of Rnits
#'
#' @name $
#' @aliases $<-,Rnits-method
#' @docType methods
#' @rdname extract-methods
#'
setReplaceMethod("$", "Rnits", function(x, name, value) {
    slot(x, name) <- value
    return(x)
}) 
