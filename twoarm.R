HC.R.BA.2A.SS.trial <- function(sdata, datah, bpar, param) {    
    sdata <- bpar$gensurv(rep.int(1, nrow(sdata)), sdata, bpar)

    ## select one historical control for each new recruit
    ##m <- max(1, nrow(datah) - nrow(sdata))
    ##datah <- datah[order(datah$entry), ][(m+1):nrow(datah),]
    
    p0 <- mean(datah$d14)
    p1 <- mean(sdata$d14)

    a <- nrow(datah)
    b <- nrow(sdata)

    s0 <- sum(datah$d14)
    s1 <- sum(sdata$d14)

    ed <- getEndpoint(rbind(datah, sdata), param$x, TRUE)

    if(ed$var[1] + ed$var[2] == 0) Z <- 0
    else {
        ##lo <- logodds(ed$n, ed$p, ed$var)
        ##if(!is.finite(lo[1])) Z <- lo[1]        
        ##else Z <- lo[1] / sqrt(lo[2])
        Z <- (ed$p[1] - ed$p[2]) / sqrt(ed$var[1] + ed$var[2])
    }
        
    su <- 1-pnorm(Z) < param$alpha
    
    c(su, max(sdata$entry + pmin(sdata$time, param$x)), ed$n[1], ed$n[2], ed$p[1], ed$p[2])
}

## parallel two-arm trials
CC.R.BA.2A.SS.par <- function(data, datah, param, bpar) {
    K <- param$K
    n <- param$n

    ## maximum available sample size per trial
    nmax <- min(n, floor(nrow(data) / K))

    ## randomly assign patients to trials
    A <- sample(rep(1:K, nmax))
    data <- data[1:length(A), ]

    res <- sapply(1:K, function(k) param$trial(data[A == k, ], datah, bpar, param))
    
    maxTime <- res[2,] + param$anaDelay
    n0 <- sum(res[3,])
    n <- res[4,]
    su <- res[1,]
    p <- res[6,]
    
    result <- any(res[1,] == 1)

    if(result) {
        index1 <- which(su == 1)
        index2 <- which.min(p[index1])
        best.trt <- (1:K)[index1][index2]
        
        dprev <- 0 ##deaths.prevented(param$x, maxTime,
                   ##               bpar$HR.D[best.trt+1], bpar$HR.R[best.trt+1],
                   ##               cumhaz=bpar$distr$cumhaz)
    } else dprev <- 0

    res <- matrix(c(rep(param$design, K+1), rep(result, K+1), 0:K, rep(1, K+1),
                    c(max(maxTime), maxTime),
                    c(n0, n), c(mean(res[5,]), res[6,]), c(FALSE, su), c(FALSE, !su),
                    rep(dprev, K+1)), nrow=K+1, ncol=10)
    
    colnames(res) <- c("design", "result", "arm", "ia", "time", "n", "p", "su", "fu", "dprev")
    res
}

## sequential two-arm trials
## CC.R.BA.2A.SS.seq <- function(data, datah, param, scm) {

##     K <- param$K
##     n <- param$n

##     K.real <- min(K, ceiling(nrow(data) / n))
##     nn <- c(rep(n, K.real - 1), min(n*K, nrow(data)) - n*(K.real - 1))

##     ## sequentially assign patients to trials (first nn[1] patients to trial 1, ...)
##     A <- rep(1:K.real, nn)
##     data2 <- data[1:length(A), ]
    
##     res <- sapply(1:K.real, function(k) param$trial(data2[A == k, ], datah, scm))

##     n0 <- sum(res[3,])
##     n <- res[4,]
    
##     maxTime <- res[2, ] + param$anaDelay

    
##     data.frame(arm=0:K, ia=rep(1, K.real),
##                time=c(max(maxTime), maxTime, rep(NA, K-K.real)),
##                n=c(n0, n, rep(NA, K-K.real)),
##                p=c(mean(res[5,]), res[6,], rep(NA, K-K.real)),
##                su=c(FALSE, res[1,], rep(FALSE, K-K.real)),               
##                fu=c(FALSE, !res[1,], rep(FALSE, K-K.real)))
## }


param2AparHC <- list(func=CC.R.BA.2A.SS.par, param=list(trial=HC.R.BA.2A.SS.trial, design=2,
                                                        n=105, dp=0.2, K=1, alpha=0.025, anaDelay=1, x=14))
