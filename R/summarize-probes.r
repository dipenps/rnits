#' Summarize probe level data to gene level data
#' 
#' The code utilizes the probe-gene mapping from features file to summarize probe-level log ratios to gene level ratios. 
#' 
#' Tukey's biweight is used to compute gene level summary
#' 
#' @export
#' @docType methods
#' @rdname summarizeProbes-methods
#' @param object \code{\linkS4class{Rnits}} object
#' @return An object of class \code{\linkS4class{Rnits}} with gene level log ratios, which can be retrieved by \code{getLR(object)}
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast  chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
#' \donttest{
#' # Summarize gene-level data
#' rnitsobj <- summarizeProbes(rnitsobj)
#'}
#'
#'
setGeneric("summarizeProbes", function(object) {
    standardGeneric("summarizeProbes")
})
#' @rdname summarizeProbes-methods
#' @aliases summarizeProbes,character,ANY-method
#' @importFrom affy tukey.biweight
setMethod("summarizeProbes", "Rnits", function(object) {
    ## Check object
    if (class(object) != "Rnits") 
        stop("Object not of Rnits class")
    ## Called Internally by fit()
    if (!"GeneName" %in% colnames(fData(object))) 
        stop("Cannot summarize probes without GeneName column.\n")
    genenames <- as.character(fData(object)[["GeneName"]])
    emptynames = which(genenames == "" | is.na(genenames))
    cat(length(emptynames), " probes had no gene names. Removing them from further processing.\n")
    if (length(emptynames) > 0) {
        newobject <- object[-emptynames, ]
    } else {
        newobject <- object
    }
    genelist <- as.character(fData(newobject)$GeneName)
    ugenes = unique(genelist)
    if (any(is.na(ugenes)) | any(ugenes == "")) 
        stop("Missing or NA values not permissible in gene names.")
    cat("Found ", length(ugenes), " unique genes. \n")
    lr = exprs(newobject)
    sinkfile = "/dev/null"
    if (.Platform$OS.type == "windows") {
        sinkfile = "NUL"
    }
    sink(sinkfile)
    if (any(is.na(lr))) 
        lr <- suppressWarnings(impute.knn(lr)$data)
    sink()
    Ng = length(ugenes)
    gene.lr = matrix(0, nrow = Ng, ncol = ncol(lr))
    gene.flag = rep(0, Ng)
    if (class(genelist) == "data.frame") 
        genelist = unlist(genelist)
    for (i in 1:Ng) {
        idx = which(genelist == ugenes[i])
        probe.dat = matrix(lr[idx, ], ncol = ncol(lr))
        gene.dat = apply(probe.dat, 2, function(x) tukey.biweight(x))
        gene.lr[i, ] = gene.dat
        gene.flag[i] = sum(fData(newobject)$Flag[idx])
    }
    row.names(gene.lr) <- ugenes
    ugeneDf <- data.frame(ID = ugenes, Flag = gene.flag)
    rownames(ugeneDf) <- ugenes
    newobject@geneData <- new("ExpressionSet", exprs = gene.lr, phenoData = phenoData(newobject), 
        featureData = AnnotatedDataFrame(data = ugeneDf))
    return(newobject)
}) 
