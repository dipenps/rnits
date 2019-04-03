# ############### All Functions
ranknormalization <- function(ch1, ch2, plot = TRUE) {
    # ch1[ch1==0] = NA ch2[ch2==0] = NA Impute with k=10
    sinkfile = "/dev/null"
    if (.Platform$OS.type == "windows") {
        sinkfile = "NUL"
    }
    sink(sinkfile)
    ch1imp <- suppressWarnings(impute.knn(as.matrix(ch1))$data)
    ch2imp <- suppressWarnings(impute.knn(as.matrix(ch2))$data)
    sink()
    # ch1imp <- ch1 ch2imp <- ch2
    logch1 = log2(ch1imp)
    logch2 = log2(ch2imp)
    ######################################################## SMOOTHING Smooth sorted average ch1 intensity and sorted ch1 intensity from
    ######################################################## individual arrays using Spline smoothing Use the fit model of the smoothing, to
    ######################################################## predict smoothed sorted ch2 intensity for individual array Median of ch1
    ######################################################## intensities across all arrays
    ch1Ref <- apply(logch1, c(1), median, na.rm = TRUE)
    ## Sort channel 1 average intensities
    refsort = order(ch1Ref)
    Pch2 <- c()
    Pch1 <- c()
    Sch1 <- c()
    ## 1. Fit smoothing function to individual array ch1 2. Adjust individual array
    ## ch1 using fit function 3.  Predict individual array ch2 using the same fit
    ## function
    for (i in 1:ncol(ch1)) {
        # print(paste('Smoothing Array', i));
        Sch1 <- smooth.spline(logch1[refsort, i], ch1Ref[refsort], df = 5)
        predRef <- predict(Sch1$fit, logch1[refsort, i])
        Pch1 <- cbind(Pch1, predRef$y)
        predSample <- predict(Sch1$fit, logch2[refsort, i])
        Pch2 <- cbind(Pch2, predSample$y)
    }
    ## Sort gene id based on intensity
    mPch1 = colMeans(Pch1)
    mPch2 = colMeans(Pch2)
    nPch = mPch2 - mPch1
    Ng = nrow(ch1)
    Nc = ncol(ch1)
    Pch2 = Pch2 - matrix(rep(nPch, Ng), Ng, Nc, byrow = TRUE)
    Pch2out = Pch2
    Pch1out = Pch1
    Pch1out[refsort, ] <- Pch1
    Pch2out[refsort, ] <- Pch2
    if (plot == TRUE) {
        par(mfrow = c(2, 1))
        boxplot(Pch1out, main = "Normalized Ch 1", cex = 0.1)
        boxplot(Pch2out, main = "Normalized Ch 2", cex = 0.1)
    }
    result = list(Pch1 = Pch1out, Pch2 = Pch2out, logratio = Pch2out - Pch1out)
    return(result)
}
solvemat <- function(Y = NULL, x = NULL) {
    Yhat <- calcYhat(Y, x)
    hatX <- hat(x)
    res <- studentize((Y - Yhat), hatX)
    return(res)
}
calcYhat <- function(Y = NULL, x = NULL) {
    coef = solve.qr(qr(x), t(Y))
    Yhat = crossprod(coef, t(x))
    return(Yhat)
}
studentize <- function(resmat, hat) {
    normf = sqrt(1 - hat)
    normf.mat = matrix(normf, nrow = nrow(resmat), ncol = ncol(resmat), byrow = TRUE)
    normres = resmat/normf.mat
}
## http://www.utulsa.edu/microarray/Articles/sam%20manual.pdf
finds0 <- function(rationum, rssALT) {
    if (any(table(rssALT) > 1)) {
        cat("Warning: Non-unique residuals observed. Check original data.\n")
        rssALT <- rssALT + rnorm(length(rssALT)) * 1e-05
    }
    s = rssALT
    qs = quantile(s, seq(0, 1, 0.01))
    alpha = seq(0, 1, 0.05)
    madval = matrix(0, 21, 99)
    quantmat = matrix(rep(quantile(s, alpha), length(s)), byrow = TRUE, nrow = length(s))
    r = apply(quantmat, 2, function(x) rationum/(s + x))
    for (k in 1:21) {
        for (i in 1:99) {
            madval[k, i] = mad(r[which(s >= qs[i] & s < qs[i + 1]), k])
        }
    }
    cvs = apply(madval, 1, function(x) {
        sd(x)/mean(x)
    })
    s0 = quantile(rssALT, alpha[which.min(cvs)])
}
getp <- function(lr, lr0) {
    m <- length(lr)
    B <- length(lr0)/m
    v <- c(rep(TRUE, m), rep(FALSE, m * B))
    v <- v[rev(order(c(lr, lr0)))]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v == TRUE] - w)/(B * m)
    p <- p[rank(-lr)]
    return(p)
}
modelfit <- function(comb = NULL, X = NULL, breaks = NULL, inp_s0 = 0) {
    Nsets = length(breaks)
    XX = c()
    for (i in 1:Nsets) {
        XX = rbind(XX, X)
    }
    hatX = hat(X)
    hatXX = hat(XX)
    # NULL calculations
    nulldata = as.matrix(comb)
    resNULL = solvemat(nulldata, XX)
    arrayid = c(breaks, ncol(comb) + 1)
    # ALT calculations
    resALT = c()
    for (i in 1:Nsets) {
        altdata = as.matrix(comb[, arrayid[i]:(arrayid[i + 1] - 1)])
        resALT = cbind(resALT, solvemat(altdata, X))
    }
    # Compute s0 and stats
    rssALT = rowSums(resALT^2)
    rssNULL = rowSums(resNULL^2)
    if (inp_s0 == 0) {
        s0 = finds0(rssNULL - rssALT, rssALT)
    } else {
        s0 = inp_s0
    }
    ratio = (rssNULL - rssALT)/(rssALT + s0)
    list(ratio = ratio, resALT = resALT, s0 = s0)
}
## Compute ratio statistic on observed and permuted datasets
tsFit <- function(inpdata = NULL, X = NULL, breaks = NULL, s0 = 0, verbatim = FALSE, 
    B = 100) {
    boot = FALSE
    Nsets = length(breaks)
    XX <- c()
    for (i in 1:Nsets) {
        XX <- rbind(XX, X)
    }
    nullHat = calcYhat(inpdata, XX)
    if (!boot) {
        fit = modelfit(inpdata, X, breaks, s0)
        ratio = fit$ratio
        resALT = fit$resALT
        s0 = fit$s0
        boot = TRUE
    }
    if (boot) {
        nboot = B
        ratiostar = matrix(0, nrow = nrow(inpdata), ncol = nboot)
        if (verbatim) 
            pb <- txtProgressBar(min = 0, max = nboot, style = 3)
        for (i in 1:nboot) {
            if (verbatim) 
                setTxtProgressBar(pb, i)
            v = sample(ncol(nullHat), replace = TRUE)
            nullcomb = nullHat + resALT[, v]
            fitb = modelfit(nullcomb, X, breaks, s0)
            ratiostar[, i] = fitb$ratio
        }
        if (verbatim) 
            close(pb)
    }
    result = list(ratio = fit$ratio, ratiostar = ratiostar, s0 = fit$s0)
    result
}
placeKnots <- function(cdata, deg, df, timevec) {
    nc = ncol(cdata)
    ng = nrow(cdata)
    nknots = df - 1 - deg
    if (nknots < 1) {
        return(NULL)
    }
    diffvec = c()
    for (k in 3:(nc - 2)) {
        diffvec[k - 2] = (mean(cdata[, k] - cdata[, k - 1]) * mean(cdata[, k] - cdata[, 
            k + 1]))
    }
    out <- timevec[order(diffvec, decreasing = TRUE)[1:nknots] + 2]
    if(any(out==max(timevec))) out[out==max(timevec)] <- max(timevec)-1
    if(any(out==min(timevec))) out[out==min(timevec)] <- min(timevec)+1 # avoid edge cases
    return(out)
} 
