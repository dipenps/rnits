#' Plot profiles of top genes/probes
#' 
#' After \code{fit} has been applied on \code{\linkS4class{Rnits}} object, plot the profiles 
#' of N top ranking genes/probes. 
#' 
#' @param object \code{\linkS4class{Rnits}} object.
#' @param id Names of specific genes or probes to be plotted. Overrides \code{fdr} and \code{top} argument.
#' @param fdr FDR cut-off plotting top probes or genes. Overrides \code{top} argument.
#' @param top Number of top genes or probes whose profile is to be plotted. Default \code{48}.
#' @param pdf Save plot as pdf? Default \code{FALSE}.
#' @param sort.by Criteria for sorting top genes or probes. Default \code{'p-value'}.
#' @param filename Name of pdf file. Default \code{topplots.pdf}.
#' @param scale_y If 'free', use free scales for plots. Default NULL.
#' @param \dots Optional arguments to plot
#' 
#'
#' @export
#' @rdname plotResults-methods
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
#' # Plot top results
#' plotResults(rnitsobj, top = 16)
#'
#'}
#'
#'
#'
setGeneric("plotResults", function(object, id = NULL, fdr = NULL, top = 48, pdf = FALSE, 
    sort.by = "p-value", filename = "TopPlots.pdf", scale_y = NULL) {
    standardGeneric("plotResults")
})
#' @rdname plotResults-methods
#' @aliases plotResults,character,ANY-method
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
setMethod("plotResults", signature = "Rnits", function(object, id = NULL, fdr = NULL, 
    top = 48, pdf = FALSE, sort.by = "p-value", filename = "TopPlots.pdf", scale_y = NULL) {
    # browser() Check object
    if (class(object) != "Rnits") 
        stop("Object not of Rnits class")
    if (!object@callData$fit.model) 
        stop("Run fit() on data first")
    if (object@callData$gene.level) {
        lr <- exprs(object@geneData)
        genes <- fData(object@geneData)$ID
    } else {
        lr <- exprs(object)
        genes <- featureNames(object)
    }
    timevec = sort(unique(pData(object)$Time))
    colnames(lr) <- paste0(pData(object)$Sample, "_", timevec)
    ## Check sort.by
    if (is.null(sort.by)) {
        sort.by = "p-value"
    }
    if (!sort.by %in% c("p-value", "FDR")) 
        stop("sort.by argument may take the value 'p-value' or 'FDR'")
    ## Check pdf
    if (is.null(pdf)) {
        pdf = FALSE
    }
    if (pdf & is.null(filename)) {
        filename = "TopPlots.pdf"
    }
    ## Check id
    if (is.null(fdr) & is.null(top) & is.null(id)) {
        top = 48
    }
    if (!is.null(id)) {
        top = NULL
        topNames <- id
        hits <- match(toupper(topNames), toupper(genes))
        if (min(!is.na(hits)) == 0) 
            stop("One or more gene/probe id's not found in data. Check id's")
        stable <- data.frame(Name = id, ClusterID = 0, Rank = 0)
    }
    if (is.null(id) & is.null(fdr) & !is.null(top)) {
        stable <- summary(object, top = top, sort.by = sort.by)
        topNames <- as.vector(stable$Name)
        hits <- match(toupper(topNames), toupper(genes))
    }
    if (is.null(id) & !is.null(fdr)) {
        stable <- summary(object, fdr = fdr, sort.by = sort.by)
        topNames <- as.vector(stable$Name)
        hits <- match(toupper(topNames), toupper(genes))
    }
    gene.level <- object@callData$gene.level
    # flag <- object@probeData$Flag if(gene.level) flag <- object@fitData$gene.flag
    ## Plot genes
    if (pdf) {
        pdf(filename)
    }
    n.win = ceiling(length(hits)/16)
    for (pnum in 1:n.win) {
        idx <- hits[((pnum - 1) * 16 + 1):min(length(hits), pnum * 16)]
        data = melt(lr[idx, ])
        data$Sample <- sapply(levels(unclass(data$Var2))[unclass(data$Var2)], function(x) strsplit(x, 
            "_")[[1]][1])
        data$Time <- as.numeric(sapply(levels(unclass(data$Var2))[unclass(data$Var2)], 
            function(x) strsplit(x, "_")[[1]][2]))
        colnames(data)[1:2] <- c("Name", "SampleTime")
        data$Cluster <- stable$ClusterID[match(data$Name, stable$Name)]
        data$Name <- paste0(data$Name, " (R", stable$Rank[match(unique(data$Name), 
            stable$Name)], ":C", stable$ClusterID[match(unique(data$Name), stable$Name)], 
            ")")
        data$Name <- factor(data$Name, levels = unique(data$Name))
        p <- ggplot(data, aes(Time, value)) + geom_point(aes(col = Sample)) + geom_smooth(aes(group = Sample, 
            col = Sample), method = "loess", alpha = 0) + ylab("Expression (log2)") + 
            theme_bw() + theme(legend.position = "top") + scale_color_brewer(type = "div", 
            palette = "Set1")
        if (!is.null(scale_y) && scale_y == "free") {
            p <- p + facet_wrap(~Name, scales = "free")
        } else {
            p <- p + facet_wrap(~Name)
        }
        suppressWarnings(print(p))
        if (!pdf) 
            reply <- readline(prompt = paste("Next:", pnum, "of", n.win, sep = " "))
    }
    if (pdf) {
        dev.off()
        cat("Plots saved to TopPlots.pdf in the working directory \n")
    }
}) 
