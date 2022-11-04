#' Input the RGlist raw data, build a Rnits object and perform filtering and normalization
#' 
#' This function takes high-dimensional expression data as a RGList, creates a 
#' \code{\linkS4class{Rnits}} object for subsequent filtering and normalization
#' 
#' See the Limma User's Guide for more details on \code{read.maimages}, \code{normalizeBetweenArrays}, 
#' \code{normalizeWithinArrays} and \code{RGList}. For importing microarray raw data, 
#' use the 'Targets file' to specify experimental design. The target file has columns 
#' SlideNumber, FileName, Cy3 (description of Cy3 channel ref/control/treatment), Cy5 
#' (description of Cy3 channel ref/control/treatment) and Time. Time values should be 
#' identical for control and treatment. 
#' 
#' @param obj Raw expression data in \code{RGlist}, \code{AffyBatch} or simple data 
#' frame format
#' @param probedata A data frame containing the probe names that should match the probe 
#' names in raw data (optional)
#' @param phenodata A data frame with information about sample names. The rownames of 
#' the data frame must match column names of the expression values. If input data is 
#' data frame of log ratios, this is required.
#' @param filter An argument to perform background filtering of probes. If \code{NULL}, 
#' no filtering is done. If an integer (0-500), probes are flagged based on raw channel 
#' intensity. If a vector of two numbers is provided, the first will be used for red channel 
#' and the second for green channel. If \code{'background'}, probes whose intensities are 
#' lower than 2 standard deviations less than the mean of the background intensity for 
#' the channel are flagged.
#' @param normalize Character string specifying the normalization method for raw data. 
#' If \code{Intensity}, the reference channels for all arrays are used to construct an 
#' array-specific smoothing function which is then applied to normalize the sample channel. 
#' If \code{Between}, the normalization method \code{normalizeBetweenArrays} in the 
#' LIMMA package is used (use \code{normmethod} to further specify normalization 
#' methods. See packaged LIMMA for details.). If \code{Within}, the normalization 
#' method \code{normalizeWithinArrays} in the LIMMA package is used.
#' @param normmethod Normalization method for input data. Default \code{NULL}. 
#' Can be one of 'quantile', 'vsn', 'Between'
#' @param background Only for AffyBatch data. If \code{TRUE}, background filtering 
#' will be done on Affy data.
#' @param center If \code{TRUE}, the log-ratio data will be mean centered to 
#' 0 in the column space. 
#' @param plot If \code{TRUE}, boxplots of normalized channel intensities and 
#' log-ratios are drawn.
#' @param threshold Default \code{0.8}. Fraction of samples with missing data 
#' for individual probes to be filtered out.
#' @param logscale Default \code{FALSE}. Is the data in logscale? If FALSE, 
#' log2 transformation is done on the data.
#' 
#' @return An object of S4 class \code{\linkS4class{Rnits}} (which is 
#' derived from class \code{\linkS4class{exprSet}}), containing the probe data, 
#' design data, expression data, phenotypical data (i.e. Time).
#' @export 
#' @docType methods
#' @seealso \code{ExpressionSet}
#' 
#' @importFrom limma normalizeBetweenArrays normalizeWithinArrays normalizeQuantiles normalizeVSN
#' @importFrom affy rma
#' @importFrom impute impute.knn
#' @import Biobase
#' @importFrom graphics boxplot
#' 
#' @examples
#' # load pre-compiled expressionSet object for Ronen and Botstein yeast chemostat data
#' data(yeastchemostat)
#' rnitsobj = build.Rnits(yeastchemostat, logscale = TRUE, normmethod = 'Between')
build.Rnits <- function(obj, probedata = NULL, phenodata = NULL, filter = NULL, normalize = NULL, 
                        normmethod = NULL, plot = FALSE, center = FALSE, background = NULL, threshold = 0.8, 
                        logscale = FALSE) {
  callData <- list()
  fdata <- pdata <- AnnotatedDataFrame(data = data.frame())
  ## Check necessary arguments
  if (!"ExpressionSet" %in% class(obj)) {
    if (is.null(phenodata)) 
      stop("Phenotype data table must be provided for object\n")
    if (!"Sample" %in% colnames(phenodata)) 
      stop("Phenotype data table must have \"Sample\" column\n")
    if (is.null(probedata)) 
      cat("Probedata not provided\n")
  }
  if (!is.null(probedata) & is.null(probedata$ProbeID)) 
    stop("Probe data table has no column named \"ProbeID")
  if ("data.frame" %in% class(obj) | "matrix" %in% class(obj)) {
    cat("Log scale is", logscale, "\n")
    ## Data frame checks
    if (nrow(phenodata) != ncol(obj)) 
      stop("Phenotype data table must have same number of columns as data table\n")
    if (all(rownames(phenodata) != colnames(obj))) 
      stop("Phenotype data table must have the same column names as data table\n")
    if (!is.null(probedata) && nrow(probedata) != nrow(obj)) 
      stop("Probe data table must equal number of rows as data table\n")
    ## Data matrix checks
    if (!logscale && min(obj, na.rm = TRUE) < 0) 
      stop("Data is indicated to be NOT in log scale, but negative values observed")
    if (!logscale) 
      obj <- log2(obj)
    callData$Type = "matrix"
    if (is.null(normmethod)) {
      cat("No filtering or normalization done for input object\n")
    } else {
      if (!normmethod %in% c("quantile", "vsn", "Between")) 
        stop("Invalid normmethod argument for data frame. Must be vsn, quantile or Between")
      if (normmethod == "quantile") {
        cat("Quantile normalization\n")
        obj <- normalizeQuantiles(obj)
      }
      if (normmethod == "vsn") {
        cat("VSN normalization\n")
        obj <- normalizeVSN(obj)
      }
      if (normmethod == "Between") {
        cat("Between Array normalization\n")
        obj <- normalizeBetweenArrays(data.matrix(obj))
      }
    }
    if (!is.null(probedata)) {
      rownames(obj) <- probedata$ProbeID
      rownames(probedata) <- probedata$ProbeID
      fdata <- AnnotatedDataFrame(data = probedata)
    } else {
      probedata <- data.frame(ProbeID = paste0("p", 1:nrow(obj)))
      fdata <- AnnotatedDataFrame(data = probedata)
      if (is.null(rownames(obj))) 
        rownames(obj) <- probedata$ProbeID
    }
    pdata <- AnnotatedDataFrame(data = phenodata)
    to <- new("Rnits", exprs = obj, phenoData = pdata, featureData = fdata, callData = callData)
    ## Probe filter
    if (!is.null(filter)) {
      if (!filter %in% c(TRUE, FALSE)) 
        stop("filter argument must be TRUE or FALSE")
      if (filter) {
        cat("Probe filter is TRUE. Applying probe threshold of ", 100 * threshold, 
            "% samples\n")
        probefilt <- which(apply(exprs(to), 1, function(x) sum(is.na(x))) > 
                             threshold * ncol(exprs(to)))
        cat(length(probefilt), " probes have NA values in more than ", ceiling(threshold * 
                                                                                 ncol(exprs(to))), "samples\n")
        if (length(probefilt) > 0) 
          to <- to[-probefilt, ]
      }
    }
  } else if (class(obj) == "RGList") {
    callData$Type = "RGList"
    ## RGlist checks
    if (nrow(phenodata) != ncol(obj)) 
      stop("Phenotype data table must have same number of columns as data table\n")
    if (!is.null(probedata) && nrow(probedata) != nrow(obj)) 
      stop("Probe data table must have equal number of rows as data table\n")
    if (!is.null(probedata) & !all(probedata$ProbeID == obj$genes$ID)) 
      cat("Probe ID's do not match RGList gene ID's\n")
    ## Filter
    if (is.null(filter)) {
      callData$Filter = "None"
      cat("No probe filtering\n")
    } else {
      if (filter[1] == "background") {
        nacount <- rep(0, nrow(obj))
        cat("Performing background filtering on data\n")
        callData$Filter = "Background"
        for (i in 1:ncol(obj$R)) {
          r.id = which(obj$R[, i] <= median(obj$Rb[, i]) + 2 * sd(obj$Rb[, 
                                                                         i]))
          g.id = which(obj$G[, i] <= median(obj$Gb[, i]) + 2 * sd(obj$Gb[, 
                                                                         i]))
          nacount[r.id] <- nacount[r.id] + 1
          nacount[g.id] <- nacount[g.id] + 1
        }
      } else if (is.numeric(filter)) {
        callData$Filter = "Intensity"
        if (max(filter) > 500 | min(filter) <= 0) 
          stop("Filter intensity should be between 0 and 500\n")
        if (length(filter) == 2) {
          filtR = filter[1]
          filtG = filter[2]
        } else if (length(filter) == 1) {
          filtR = filtG = filter
        } else stop("filter intensity should be either one or two integers")
        nacount <- rep(0, nrow(obj))
        cat("Performing Intensity filtering on data\n")
        for (i in 1:ncol(obj$R)) {
          r.id = which(obj$R[, i] <= filtR)
          g.id = which(obj$G[, i] <= filtG)
          nacount[r.id] <- nacount[r.id] + 1
          nacount[g.id] <- nacount[g.id] + 1
        }
      } else stop("Filter argument should either be 'background' or an intensity value between 0 and 500")
      # Probe filter
      thr = threshold
      thr = ncol(obj$R) * thr
      filt = which(nacount <= thr)
      cat(paste0("\tProbes before filtering: ", nrow(obj)), "\n")
      obj = obj[filt, ]
      if (!is.null(probedata)) 
        probedata = probedata[filt, ]
      cat(paste0("\tProbes after filtering: ", nrow(obj), "\n"))
    }
    ## Normalize
    if (is.null(normalize)) {
      callData$Normalize = "None"
      cat("No normalization done\n")
      R <- obj$R
      G <- obj$G
      lr <- log2(R/G)
    } else {
      if (normalize == "Between") {
        callData$Normalize = "Between"
        cat("Between Array Normalization from Limma\n")
        if (is.null(normmethod)) {
          normmethod = "Aquantile"
        } else if (!normmethod %in% c("none", "scale", "quantile", "Aquantile", 
                                      "Gquantile", "Rquantile", "Tquantile", "vsn")) {
          cat("Normalization Method ", normmethod, " not a valid option for Between Array normalization. Continuing with Aquantile\n")
          normmethod = "Aquantile"
        }
        callData$NormMethod = normmethod
        MAobj <- normalizeBetweenArrays(obj, method = normmethod)
        lr <- MAobj$M
        R <- MAobj$A + lr/2
        G <- R - lr
      } else if (normalize == "Within") {
        callData$Normalize = "Within"
        cat("Within Normalization from Limma\n")
        if (is.null(normmethod)) {
          normmethod = "loess"
        } else if (!normmethod %in% c("none", "median", "loess", "printtiploess", 
                                      "composite", "control", "robustspline")) {
          cat("Normalization Method ", normmethod, " not a valid option for Within Array normalization. Continuing with loess.\n")
          normmethod = "loess"
        }
        callData$NormMethod = normmethod
        MAobj <- normalizeWithinArrays(obj, method = normmethod)
        lr <- MAobj$M
        R <- MAobj$A + lr/2
        G <- R - lr
      } else if (normalize == "Intensity" | normalize == "IntensityGreen") {
        cat("Intensity Normalization on Green Channel\n")
        callData$Normalize = "Intensity"
        outList <- ranknormalization(ch1 = obj$G, ch2 = obj$R, plot = FALSE)
        R <- outList$Pch2
        G <- outList$Pch1
        lr <- outList$logratio
      } else if (normalize == "IntensityRed") {
        cat("Intensity Normalization on Red Channel\n")
        outList <- ranknormalization(ch1 = obj$R, ch2 = obj$G, plot = FALSE)
        R <- outList$Pch1
        G <- outList$Pch2
        lr <- -outList$logratio
      } else stop("Invalid normalization argument")
    }
    obj$R <- R
    obj$G <- G
    if (!is.null(probedata)) {
      rownames(lr) <- rownames(obj$R) <- rownames(obj$G) <- probedata$ProbeID
      rownames(probedata) <- probedata$ProbeID
      fdata <- AnnotatedDataFrame(data = probedata)
    } else {
      probedata <- obj$genes
      rownames(lr) <- rownames(obj$R) <- rownames(obj$G) <- rownames(probedata) <- obj$genes$ID
      fdata <- AnnotatedDataFrame(data = probedata)
    }
    pdata <- AnnotatedDataFrame(data = phenodata)
    to <- new("Rnits", rawTwoColor = obj, exprs = lr, phenoData = pdata, featureData = fdata, 
              callData = callData)
  } else if (class(obj) == "AffyBatch") {
    calldata$Type = "AffyBatch"
    nmOpt <- bgOpt <- TRUE
    if (is.null(normalize)) {
      cat("Setting normalization ON for Affy object. Disable by settingn normalize = FALSE\n")
      nmOpt <- TRUE
    }
    if (is.null(background)) {
      cat("Setting background correction ON for Affy object. Disable by setting background = FALSE\n")
      bgOpt <- TRUE
    }
    rmaset <- rma(obj, verbose = FALSE, background = bgOpt, normalize = nmOpt)
    to <- new("Rnits", rawAffy = obj, exprs = exprs(rmaset), featureData = featureData(rmaset), 
              phenoData = phenoData(rmaset), callData = callData)
  } else if (class(obj) == "ExpressionSet") {
    callData$Type = "ExpressionSet"
    phenodata = pData(obj)
    featuredata = fData(obj)
    if (!"Sample" %in% colnames(phenodata)) 
      stop("Phenotype data table must have Sample column\n")
    if (!"Time" %in% colnames(phenodata)) 
      stop("Phenotype data table must have Time column\n")
    dat <- exprs(obj)
    if (!logscale && min(dat, na.rm = TRUE) < 0) 
      stop("Data is indicated to be NOT in log scale, but negative values observed")
    if (!logscale) 
      dat <- log2(dat)
    callData$Type = "ExpressionSet"
    if (is.null(normmethod)) {
      cat("No filtering or normalization done for input object\n")
      exprs(obj) <- dat
    } else {
      if (!normmethod %in% c("quantile", "vsn", "Between")) 
        stop("Invalid normmethod argument for data frame. Must be vsn, quantile or Between")
      if (normmethod == "quantile") {
        cat("Quantile normalization\n")
        dat <- normalizeQuantiles(dat)
      }
      if (normmethod == "vsn") {
        cat("VSN normalization\n")
        dat <- normalizeVSN(dat)
      }
      if (normmethod == "Between") {
        cat("Between Array normalization\n")
        dat <- normalizeBetweenArrays(data.matrix(dat))
      }
      exprs(obj) <- dat
    }
    if (!"GeneName" %in% colnames(featuredata)) {
      if (!"Gene Symbol" %in% colnames(featuredata)) {
        stop("GeneName or Gene Symbol column not found in Probe data.\n")
      } else {
        featuredata[["GeneName"]] <- as.character(featuredata[["Gene Symbol"]])
      }
    }
    to <- new("Rnits", exprs = exprs(obj), featureData = as(featuredata, "AnnotatedDataFrame"), 
              phenoData = as(phenodata, "AnnotatedDataFrame"), callData = callData)
    ## Probe filter
    if (!is.null(filter)) {
      if (!filter %in% c(TRUE, FALSE)) 
        stop("filter argument must be TRUE or FALSE")
      if (filter) {
        cat("Probe filter is TRUE. Applying probe threshold of ", 100 * threshold, 
            "% samples\n")
        probefilt <- which(apply(exprs(to), 1, function(x) sum(is.na(x))) > 
                             threshold * ncol(exprs(to)))
        cat(length(probefilt), " probes have NA values in more than ", ceiling(threshold * 
                                                                                 ncol(exprs(to))), "samples\n")
        if (length(probefilt) > 0) 
          to <- to[-probefilt, ]
      }
    }
  } else stop("Input obj must be either a data frame, matrix, RGList, ExpressionSet or AffyBatch.")
  
  pdata <- pData(to)
  
  ## Replicate sample data
  samp_time <- with(pdata, paste0(Sample, "_", Time))
  if (max(table(samp_time)) > 1) {
    #stop("Some series have replicate samples.")
    #print("Some series have replicate samples.")
    if(length(unique(table(samp_time))) > 1) stop("Uneven number of replicates")
    if(!'Replicate' %in% colnames(pdata)) stop("Phenotype data must have 'Replicate' column for replicate data")
    Nrep=unique(table(samp_time))
    sort_order <- order(pdata$Sample,  pdata$Replicate, as.numeric(pdata$Time))
    rep=TRUE
  } else{
    Nrep=1
    rep=FALSE
    sort_order <- order(pdata$Sample, as.numeric(pdata$Time))
  }
  
  ## Sort sample data
  to <- to[, sort_order]
  pdata <- pData(to)
  
  
  utime <- unique(as.numeric(pdata$Time))
  usamp <- sort(unique(pdata$Sample))
  Nsets <- length(usamp) - 1
  Nc <- length(utime)*Nrep
  breaks = 1
  for (i in 1:Nsets) breaks[i + 1] = breaks[i] + Nc
  
  if ("control" %in% tolower(usamp)) {
    control_samp <- which(tolower(pdata$Sample) %in% "control")
  } else {
    control_samp <- 1:Nc
  }
  
  callData$Nsets <- Nsets
  callData$Nc <- Nc
  callData$breaks <- breaks
  callData$sample_names <- usamp
  callData$time_vec <- sort(utime)
  callData$control_samp <- control_samp
  callData$Nrep <- Nrep
  to$callData <- callData
  cat("Dataset with ", nrow(to), " features and ", Nsets + 1, " samples in ", Nrep, "replicates\n")
  
  if (center) {
    cat("Centering data\n")
    globalmean <- mean(colMeans(exprs(to), na.rm = TRUE))
    exprs(to) <- globalmean + scale(exprs(to), center = TRUE, scale = FALSE)
  }
  if (plot) 
    boxplot(exprs(to), main = "Expression", xlab = "Samples", cex = 0.5, pch = 19)
  return(to)
} 
