#' Get Normalized channel data for two channel arrays
#' 
#' For two color data, extract normalized channel data from \code{\linkS4class{Rnits}} object
#' 
#' @param object \code{\linkS4class{Rnits}} object
#' @export
#' @docType methods
#' @rdname getNormTwoChannel-methods
#' @return A list containing R and G fields for normalized Red and Green channel data respectively.

setGeneric("getNormTwoChannel", function(object) {
    standardGeneric("getNormTwoChannel")
})
#' @rdname getNormTwoChannel-methods
#' @aliases getNormTwoChannel,character,ANY-method
setMethod("getNormTwoChannel", signature = "Rnits", function(object) {
    if (object@callData$Type != "RGList") {
        stop("Data must be an RGList to retrieve normalized data.")
    }
    return(list(R = object@rawTwoColor$R, G = object@rawTwoColor$G))
}) 
