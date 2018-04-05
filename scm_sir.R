library(ggplot2)

source("util.R")
source("dcm_sir.R")

loadCxx <- function() {
    require(Rcpp)
    Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
    sourceCpp("scm_seird.cpp", cacheDir="/tmp/rcpp")
}

scmSEIRD <- function(param) {
    x <- scmSEIRDcxx(c(param$beta1, param$beta2), c(param$delta1, param$delta2),
                     param$gamma, c(param$f1, param$f2), param$changePoint,
                     param$tbegin, param$tend, param$N, param$E0, param$I0)
    ##x <- round(x, 0)
    colnames(x) <- c("time", "S", "E", "I", "R", "D", "newI")
    x <- as.data.frame(x)
    class(x) <- c("scmSEIRD", "data.frame")
    x
}

prob <- function(x, t1, par) {
    dcm <- SEIRD(par)
    
    CP <- par$changePoint
    f1 <- par$f1 * par$HR.D1[2]
    f2 <- par$f2 * par$HR.D2[2]
    beta1 <- par$beta1
    beta2 <- par$beta2
    delta1 <- par$delta1
    delta2 <- par$delta2
    gamma <- par$gamma
    
    F <- function(t, s) f1*(1 - exp(-delta1*s))*(s < (CP - t)) + (f1*(1 - exp(-delta1*max(0, CP-t))) + f2*(exp(-delta2*max(0, CP-t)) - exp(-delta2*s)))*(s >= (CP - t))

    h1 <- beta1*dcm$S / dcm$S[1] * dcm$I + gamma*dcm$E + delta1*dcm$I
    h2 <- beta2*dcm$S / dcm$S[1] * dcm$I + gamma*dcm$E + delta2*dcm$I
    h <- h1*(dcm$time < CP) + h2*(dcm$time >= CP)

    dt <- c(0, diff(dcm$time))
    g <- gamma*dcm$E*exp(-cumsum(h*dt)/dcm$S[1]) * dt / dcm$S[1]

    ## TODO: conditionally on is.finite(t1)
    sum(F(dcm$time, x) * g) / sum(g)
}


testGenScmSEIRD <- function(par) {
    x <- scm(par)

    if(is.finite(x$t1)) {
        n <- min(nrow(x$df1), 200)
        data <- genScmSEIRD(sample(c(0,1), size=n, replace=TRUE),
                            x$data1$entry[1:n], par)
        c(mean(data$d14[data$R == 0]), mean(data$d14[data$R == 1]))
    } else c(NA, NA)
}

probScm <- function(y, t1, par) {
    x <- scmSEIRD(par)

    CP <- par$changePoint
    f1 <- par$f1 * par$HR.D1[2]
    f2 <- par$f2 * par$HR.D2[2]
    beta1 <- par$beta1
    beta2 <- par$beta2
    delta1 <- par$delta1
    delta2 <- par$delta2
    gamma <- par$gamma
    
    F <- function(t, s) f1*(1 - exp(-delta1*s))*(s < (CP - t)) + (f1*(1 - exp(-delta1*max(0, CP-t))) + f2*(exp(-delta2*max(0, CP-t)) - exp(-delta2*s)))*(s >= (CP - t))

    h1 <- beta1*x$S / x$S[1] * x$I + gamma*x$E + par$delta1*x$I
    h2 <- beta2*x$S / x$S[1] * x$I + gamma*x$E + par$delta2*x$I
    h <- h1*(x$time < CP) + h2*(x$time >= CP)

    dt <- c(0, diff(x$time))
    g <- gamma*x$E*exp(-cumsum(h*dt/x$S[1])) * dt / x$S[1]

    x <- seq(0, 700, length.out=1000)
    fx <- sapply(x, function(xx) F(xx, y))

    data.frame(x=x, y=fx)
    ##print(F(t1, y))

    ##print(c(sum(g), max(cumsum(x$newI))/x$S[1]))

    ##sum(F(x$time[x$time >= t1 & x$time < (CP-y)], y) * g[x$time >= t1 & x$time < (CP-y)]) / sum(g[x$time >= t1 & x$time < (CP-y)])
}

summary.scmSEIRD <- function(x, ...) {
    z <- cumsum(x[,7])
    totalI <- z[length(z)]
    totalD <- x[nrow(x), 6]
    c(nrow(x), max(x[,1]), totalI, totalD)
}

simScmSEIRD <- function(M, param, grid, seed=NULL) {
    if(!is.null(seed)) set.seed(seed)

    ##res <- matrix(0, nrow=length(grid), ncol=4)
    
    f <- function(i) {
        x <- scmSEIRD(param)
        summary(x)
        ##fI <- approxfun(x[,1], x[,4], method="linear", rule=2, yleft=0)
        ##fD <- approxfun(x[,1], x[,6], method="linear", rule=2, yleft=0)
        ##fnewI <- approxfun(x[,1], x[,7], method="linear", rule=2, yleft=0)
        ##matrix(c(grid, fI(grid), fD(grid), fnewI(grid)), nrow=length(grid)) / M
    }

    ##for(i in 1:M) res <- res + f(i)
    replicate(M, f())
    
    ##res <- as.data.frame(res)
    ##colnames(res) <- c("time", "I", "D", "newI")
    ##class(res) <- c("simScmSEIRD", "data.frame")
    ##res
}

histScmSEIRD <- function(scm, dcm) {
    df <- data.frame(totalI=scm[3,], totalD=scm[4,])
    
    p1 <- ggplot(df, aes(totalI)) + geom_histogram(col="black") + ##, fill="grey", alpha = .75) +
        labs(x="Total Number of Infections", y="Count") + theme_bw() +
        geom_vline(aes(xintercept=max(cumsum(dcm$newI)), col=I("green"))) +
        geom_vline(aes(xintercept=mean(totalI), col=I("red")))
        ##scale_color_manual(name = "Mean", values = c(DCM = I("green"), SCM = I("red")))

    p2 <- ggplot(df, aes(totalD)) + geom_histogram(col="black") +
        labs(x="Total Number of Deaths", y="Count") + theme_bw() +
        geom_vline(aes(xintercept=max(dcm$D), col=I("green"))) +
        geom_vline(aes(xintercept=mean(totalD), col=I("red")))
        ##scale_color_manual(name = "Mean", values = c(DCM = I("green"), SCM = I("red")))

    p1 <- p1 + theme(axis.text=element_text(size=13),
                     axis.title=element_text(size=14),
                     legend.text=element_text(size=13),
                     legend.title=element_text(size=14))

    p2 <- p2 + theme(axis.text=element_text(size=13),
                     axis.title=element_text(size=14),
                     legend.text=element_text(size=13),
                     legend.title=element_text(size=14))
    
    multiplot(p1, p2, cols=2)

    list(p1, p2)
    
    ##par(mfrow=c(1, 2))
    ## hist(scm[3,], main="Histogram", xlab="Total number of infections")
    ## abline(v=mean(scm[3,]), col="blue")
    ## abline(v=max(cumsum(dcm$newI)), col="red")
    ## legend("topright", c("DCM", "Mean of SCM"), col=c("red", "blue"), lty=c(1, 1))

    ## hist(scm[4,], main="Histogram", xlab="Total number of deaths")
    ## abline(v=mean(scm[4,]), col="blue")
    ## abline(v=max(dcm$D), col="red")
    ## legend("topright", c("DCM", "Mean of SCM"), col=c("red", "blue"), lty=c(1, 1))
}

## R = vector of trt assignments (0=control, ...)
## entry = vector of entry times
## param = dcm parameter list
genScmSEIRD <- function(R, entry, param) {
    hD1 <- param$f1 * param$delta1 * param$HR.D1[R+1]
    hD2 <- param$f2 * param$delta2 * param$HR.D2[R+1]
    hR1 <- (1-param$f1 * param$HR.D1[R+1]) * param$delta1
    hR2 <- (1-param$f2 * param$HR.D2[R+1]) * param$delta2 
    
    n <- length(entry)
    U <- 1-runif(n)
    cp <- param$changePoint
    ts <- cp - entry
    ts[ts < 0] <- 0
    us <- exp(-hD1*ts)
    TD <- (U < us)*(-log(U) + (hD2 - hD1)*ts) / hD2 + (U >= us)*(-log(U) / hD1)
    
    V <- 1-runif(n)
    vs <- exp(-hR1*ts)
    TR <- (V < vs)*(-log(V) + (hR2 - hR1)*ts) / hR2 + (V >= vs)*(-log(V) / hR1)
    
    time <- pmin(TD, TR)
    status <- (TD > TR) + 2*(TD <= TR)
    d14 <- status & TD <= 14  ## death within 14 days
    d28 <- status & TD <= 28  ## death within 28 days

    data.frame(id=1:n, R=R, entry=entry, time=time, status=status,
               from=0, to=status, d14=d14, d28=d28)
}

scm <- function(param) {
    x <- param$func(param)

    cnI <- cumsum(x$newI)
    
    if(max(cnI) < param$trigger) { ## trigger threshold never passed
        t0 <- Inf
        t1 <- Inf
        df0 <- data.frame()
        df1 <- data.frame()
        dfh <- data.frame()
    } else {        
        ## t0 = time clinical trial program is triggered
        t0 <- min((1:length(cnI))[cnI >= param$trigger]) + 1

        ## t1 = recruitment start = trigger time + setup time
        t1 <- t0 + param$setup

        ## program setup takes too long
        if(t1 > max(x$time[x$time >= t0 & x$newI > 0])) {
            t1 <- Inf
            df0 <- data.frame()
            df1 <- data.frame()
            dfh <- data.frame()
        } else {
            ## only recruit newly infected (rationale: older infected already under treatment)    
            entry <- x$time[x$time >= t1 & x$newI > 0] ## entry time == infection time

            df0 <- data.frame(entry=entry)
            
            ## sample up to e.max patients on each day
            df1 <- slice(group_by(df0, entry), 1:min(param$e.max, length(entry)))

            ## historical controls (everyone treated before trial program starts)
            hentry <- x$time[x$time < t1 & x$newI > 0]
            dfh <- param$gensurv(0, hentry, param)
        }
    }
    
    res <- list(t0=t0, t1=t1, data0=df0, data1=df1, datah=dfh) ## nIt1=cnI[t1], nDt1=x$D[t1])
    class(res) <- c("scm", "list")
    res
}

summary.scm <- function(x, ...) {
    ## t0 = trigger time
    ## t1 = start of recruitment
    ## n1 = number of newly infected patients after t1
    ## n2 = number of patients that can be recruited
    data.frame(t0=x$t0, t1=x$t1, n1=nrow(x$data0), n2=nrow(x$data1), ##nIt1=x$nIt1, nDt1=x$nDt1,
               I=max(cumsum(x$scm$newI)), D=max(x$scm$D))
}
