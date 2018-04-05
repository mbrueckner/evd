library(flexsurv)
library(etm)
library(mvna)

source("data.R") ## for getTimeCourse

## evd.data <- get(load("../data/evd.RData"))

evd.cif <- function(data) {
    tra <- matrix(FALSE, nrow=3, ncol=3)
    tra[1,2:3] <- TRUE
    res <- etm(data, state.names=c("0", "1", "2"), tra=tra, cens.name="cens", s=0, covariance=FALSE)   
    data.frame(time=res$time, alive=res$est[1,2,], dead=res$est[1,3,])
}

evd.cumhaz <- function(data) {
    tra <- matrix(FALSE, nrow=3, ncol=3)
    tra[1,2:3] <- TRUE
    x <- mvna(data, state.names=c("0", "1", "2"), tra=tra, cens.name="cens")
    data.frame(time=x$time, alive=x[[1]]$na, dead=x[[2]]$na)
}

## df=output of evd.cumhaz
evd.true.cif <- function(df, HR.D=1, HR.R=1, x=NULL) {
    df <- data.frame(time=df$time,
                     p=cumsum(diff(c(0, df$dead))*HR.D*exp(-df$alive*HR.R-df$dead*HR.D)))

    if(!is.null(x)) {
        f <- approxfun(df$time, df$p, method="constant", yleft=0, rule=2, f=0)
        f(x)
    } else df
}

deaths.prevented <- function(x, cal.time, HR.D, HR.R, data=evd.data,
                             cumhaz=fsCumhaz(fspline, seq(0, 260, length.out=1000))) {
    
    ## all patients hospitalized before cal.time, treated with standard care
    data0 <- data[data$entry < cal.time, ]
    
    ## all patients hospitalized after cal.time, treated with experimental treatment
    data1 <- data[data$entry >= cal.time, ]
    
    ## mortality rate under standard care
    cf <- evd.true.cif(cumhaz, 1, 1)
    p00 <- max(cf$p[cf$time <= x])

    ## mortality rate of patients recruited before trial conclusion under standard care   
    p01 <- p00 ##max(evd.true.cif(df0, 1, 1)$p)

    ## mortality rate of patients recruited after trial conclusion under experimental trt.
    cf1 <- evd.true.cif(cumhaz, HR.D, HR.R)
    p11 <- max(cf1$p[cf1$time <= x])

    ## number of deaths everyone standard care - number of deaths standard care until trial
    ## conclusion - number of deaths experimental treatment after trial conclusion
    nrow(data)*p00 - (nrow(data0)*p00 + nrow(data1)*p11)
}

LL <- function(k, theta) -sum(dgamma(evd.data$entry, shape=k, scale=theta, log=TRUE))

    
nonpar <- function(data) {
    df <- evd.cumhaz(data)

    haz1 <- diff(c(0, df$alive))
    haz2 <- diff(c(0, df$dead))

    f1 <- approxfun(df$time, haz1, method="constant", yleft=0, rule=2, f=0)
    f2 <- approxfun(df$time, haz2, method="constant", yleft=0, rule=2, f=0)

    list(f1=f1, f2=f2, cumhaz=df)
}

flexSpline <- function(data) {    
    fs1 <- flexsurvspline(Surv(time, status==1) ~ 1, data=data, k=4)
    fs2 <- flexsurvspline(Surv(time, status==2) ~ 1, data=data, k=4)

    df <- fsCumhaz(list(fs1, fs2), seq(0, max(data$time), length.out=1000))
    
    haz1 <- diff(c(0, df$alive))
    haz2 <- diff(c(0, df$dead))

    f1 <- approxfun(df$time, haz1, method="constant", rule=2, yleft=0, f=0)
    f2 <- approxfun(df$time, haz2, method="constant", rule=2, yleft=0, f=0)

    ##df1 <- summary(fs1, type="haz")[[1]]
    ##df2 <- summary(fs2, type="haz")[[1]]

    list(fs1=fs1, fs2=fs2, cumhaz=df, f1=f1, f2=f2)
}

fsCumhaz <- function(lst, times) {
    fs1 <- lst[[1]]
    fs2 <- lst[[2]]
    
    Haz1 <- Hsurvspline(times, fs1$opt$par, knots=fs1$knots)
    Haz2 <- Hsurvspline(times, fs2$opt$par, knots=fs2$knots)

    data.frame(time=times, alive=Haz1, dead=Haz2)
}

true.cif.fs <- function(data, HR.D=1, HR.R=1, times=seq(0, 50, length.out=100), x=14) {
    lst <- flexSpline(data)
    df <- fsCumhaz(lst, times)
    evd.true.cif(df, HR.D, HR.R, x=x)
}

genSurv <- function(R, data, param) {
    stopifnot(length(R) == nrow(data))

    n <- length(R)
    chaz <- param$distr$cumhaz
    f1 <- param$distr$f1
    f2 <- param$distr$f2
    taper <- param$taper
    
    time <- chaz$time
    status <- logical(n)
    entry <- data$entry

    U <- runif(n)
    tp <- taper(R, entry)
    mt <- max(time)
    
    ##sample time-to-first-transition
    tmp <- sapply(1:n, function(i) {
         s <- -chaz$alive*tp[i,1]-chaz$dead*tp[i,2] <= log(1-U[i])
         if(any(s)) c(time[s][1], TRUE)
         else c(mt, FALSE)
    })
    
    time <- tmp[1,]
    status <- as.logical(tmp[2,])
    
    ##for(r in unique(R)) {
    ##   km <- sampleFromKM(sum(R == r), chaz$time, exp(-chaz$alive*HR.R[r+1]-chaz$dead*HR.D[r+1]))
    ##   time[R == r] <- km$time
    ##   status[R == r] <- km$status
    ##}

    evtime <- time[status]   
    evtp <- tp[(1:n)[status],]
        
    rate <- f2(evtime)*evtp[2] / (f1(evtime)*evtp[1] + f2(evtime)*evtp[2])

    fstatus <- factor(levels=c(1,2,"cens"))
    fstatus[status] <- 1 + rbinom(length(evtime), 1, rate)
    fstatus[!status] <- "cens"
    
    d14 <- (fstatus == 2) & (time <= 14)  ## death within 14 days
    d28 <- (fstatus == 2) & (time <= 28)  ## death within 28 days

    df <- data.frame(id=1:n, R=R, entry=entry, time=time, status=fstatus,
                     from=0, to=fstatus, d14=d14, d28=d28,
                     centre=data$centre, country=data$country)
    rownames(df) <- NULL
    df
}

survEvd <- function(param) {
    x <- evd.data

    df0 <- data.frame()
    df1 <- data.frame()
    dfh <- data.frame()
    t0 <- Inf
    t1 <- Inf

    ## FIXME: should be part of bpar$distr
    ## x$entry <- round(rgamma(nrow(x), shape=11, scale=19))

    if(param$trigger <= length(x$entry)) {
        t0 <- sort(x$entry)[param$trigger]
        t1 <- t0 + param$setup

        if(t1 < max(x$entry)) {
            xgt1 <- x[x$entry >= t1,]
            entry <- xgt1$entry
            
            df0 <- data.frame(entry=entry, country=xgt1$country, centre=xgt1$centre)
            
            ## sample up to e.max patients on each day
            df1 <- slice(group_by(df0, entry), 1:min(param$e.max, length(entry)))

            ## historical controls (everyone treated before trial program starts)
            xlt1 <- x[x$entry < t1,]
            hentry <- xlt1$entry
            dfh <- genSurv(rep.int(0, length(hentry)), xlt1, param)
        }
    }

    res <- list(t0=t0, t1=t1, data0=df0, data1=df1, datah=dfh)
    class(res) <- c("survEvd", "list")
    res
}

summary.survEvd <- function(x, ...) {
    ## t0 = trigger time
    ## t1 = start of recruitment
    ## n0 = number of newly infected patients after t1
    ## n1 = number of patients that can be recruited
    data.frame(t0=x$t0, t1=x$t1, n0=nrow(x$data0), n1=nrow(x$data1))
}

plotResultTimeCourse <- function(df, design=unique(df$df2$design),
                                 tc=get(load("timecourse.res"))) {
    ##tc <- getTimeCourse(data)
    ##plot(tc$time, tc$total, type="l", xlab="Time [days]", ylab="Total number conf. + hosp.")
    p <- ggplot(tc, aes(x=time, y=total)) + geom_line() + xlab("Time [days]") + ylab("Total number confirmed and hospitalized") + theme_bw()
    
    ##abline(v=df$df0$mean.t0, col="green")
    p <- p + geom_vline(aes(xintercept=df$df0$mean.t0), linetype=2)
    
    ##abline(v=df$df0$min.t0, lty=2)
    ##abline(v=df$df0$max.t0, lty=2)
    
    ##abline(v=df$df0$mean.t1, col="blue")
    p <- p + geom_vline(aes(xintercept=df$df0$mean.t1), linetype=2)
    
    ##abline(v=df$df0$min.t1, lty=2, col="blue")
    ##abline(v=df$df0$mmax.t1, lty=2, col="blue")

    df1 <- df$df1[c(1, 2, 6, 3, 4, 5),]
    df1 <- df1[df$df1$design %in% design,]
    ##abline(v=df1$mean.time, col="red", lty=1:4)

    p <- p + geom_vline(aes(xintercept=mean.time, col=design), df1)

    ##abline(v=df1$min.time, lty=2, col="red")
    ##abline(v=df1$max.time, lty=2, col="red")
    
    ##iat <- df$df3 %>% summarise(time=max(time))
    ##abline(v=iat$time[iat$design %in% design], col="blue")

    ## legend(x="bottomright",
    ##        legend=c("Trigger", "Recruitment start", "FA (Single-stage)",
    ##                 "FA (GSD)", "FA (GSD Bayes)", "FA (GSD RAR)"),
    ##        col=c("green", "blue", "red", "red", "red", "red"),
    ##        lty=c(1, 1, 1, 2, 3, 4))

    p + scale_color_hue(labels=c("TACC", "TAHC", "TA.GSD", "TA.GSD (Bayes)", "TA.RAR (Bayes)", "TACC (Bayes)")) +
        guides(color=guide_legend("Designs"))
}


plotResultTimeCourseMAMS <- function(df, design=unique(df$df2$design),
                                     tc=get(load("timecourse.res"))) {
    ##tc <- getTimeCourse(data)
    ##plot(tc$time, tc$total, type="l", xlab="Time [days]", ylab="Total number conf. + hosp.")
    p <- ggplot(tc, aes(x=time, y=total)) + geom_line() + xlab("Time [days]") + ylab("Total number confirmed and hospitalized") + theme_bw()
    
    ##abline(v=df$df0$mean.t0, col="green")
    p <- p + geom_vline(aes(xintercept=df$df0$mean.t0), linetype=2)
    ##abline(v=df$df0$min.t0, lty=2)
    ##abline(v=df$df0$max.t0, lty=2)
    
    ##abline(v=df$df0$mean.t1, col="blue")
    p <- p + geom_vline(aes(xintercept=df$df0$mean.t1), linetype=2)
    ##abline(v=df$df0$min.t1, lty=2, col="blue")
    ##abline(v=df$df0$mmax.t1, lty=2, col="blue")

    df1 <- df$df1[c(3, 1, 2), ]
    
    df1 <- df1[df$df1$design %in% design,]
    ##abline(v=df1$mean.time, col="red", lty=1:3)
    p <- p + geom_vline(aes(xintercept=mean.time, col=design), df1)
    
    ##abline(v=df1$min.time, lty=2, col="red")
    ##abline(v=df1$max.time, lty=2, col="red")
    
    ##iat <- df$df3 %>% summarise(time=max(time))
    ##abline(v=iat$time[iat$design %in% design], col="blue")

    ##legend(x="bottomright",
    ##       legend=c("Trigger", "Recruitment start", "FA (MAMS)",
    ##                "FA (MAMS Bayes)", "FA (MAMS RAR)"),
    ##       col=c("green", "blue", "red", "red", "red"),
    ##       lty=c(1, 1, 1, 2, 3))

    p + scale_color_hue(labels=c("MAMS (Bayes)", "MAMS.RAR (Bayes)", "MAMS")) +
        guides(color=guide_legend("Designs"))
}
