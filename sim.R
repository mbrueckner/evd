library(dplyr)

source("util.R")
source("evd.R")

##source("scm_sir.R")

source("driver.R")
source("rar.R")
source("twoarm.R")
source("gsd.R")
source("mams.R")

enableDebug <- function() {
    options(error=recover)
    options(warn=2)
}

disableDebug <- function() {
    options(error=NULL)
    options(warn=1)
}

if(!exists("fspline")) {
    evd.data0 <- get(load("../data/evd.RData"))
    
    evd.data <- evd.data0
    evd.np <- NA ##nonpar(evd.data)
    fspline <- flexSpline(evd.data)
}

setup <- 42

## p_0, p_1 fixed
param2AH0a.FS <- list(name="2AH0a.FS", backend=survEvd, ## country="Guinea", centre="Conakry 2"
                      bpar=list(distr=fspline,
                                taper=function(r, e) matrix(1, nrow=length(e), ncol=2),
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL, paramBayesGSD, paramBGSDRAR))

## p_0, p_1 decreasing
param2AH0b.FS <- list(name="2AH0b.FS", backend=survEvd,
                      bpar=list(distr=fspline,
                                taper=function(r, e) {
                                    l <- pmin(e/250, 1)
                                    matrix(exp(c(l*log(0.68), l*log(0.68))), ncol=2)
                                },
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                   paramBayesGSD, paramBGSDRAR))


## fixed p_0, p_1
param2AH1a.FS <- list(name="2AH1a.FS", backend=survEvd,
                      bpar=list(distr=fspline,
                                taper=function(r, e) matrix(exp(rep.int(0, 2*length(e))*(r==0)
                                                            + rep(log(0.42), 2*length(e))*(r==1)),
                                                            nrow=length(e), ncol=2),
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                   paramBayesGSD, paramBGSDRAR))

## p_0 decreasing, fixed p_1
param2AH1b.FS <- list(name="2AH1b.FS", backend=survEvd,
                      bpar=list(distr=fspline,
                                taper=function(r, e) {
                                    l <- pmin(e/250, 1)
                                    matrix(exp(c(l*log(0.68), l*log(0.68))*(r == 0) +
                                           rep(log(0.42), 2*length(e))*(r == 1)), nrow=length(e), ncol=2)
                                },
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                   paramBayesGSD, paramBGSDRAR))

## p_0, p_1 decreasing
param2AH1c.FS <- list(name="2AH1c.FS", backend=survEvd,
                      bpar=list(distr=fspline,
                                taper=function(r, e) {
                                    l <- pmin(e/250, 1)
                                    matrix(exp(c(l*log(0.68), l*log(0.68))*(r == 0) +
                                           c((1-l)*log(0.42) + l*log(0.2), (1-l)*log(0.42) + l*log(0.2))*(r == 1)),
                                           nrow=length(e), ncol=2)
                                },
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                   paramBayesGSD, paramBGSDRAR))

paramMAH0a.FS <- list(name="MAH0a.FS", backend=survEvd,
                     bpar=list(distr=fspline,
                               taper=function(r, e) matrix(1, nrow=length(e), ncol=2),
                               gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                     designs=list(paramMAMS, paramCPBAL, paramCPRAR))

paramMAH0b.FS <- list(name="MAH0b.FS", backend=survEvd,
                      bpar=list(distr=fspline,
                                taper=function(r, e) {
                                    l <- pmin(e/250, 1)
                                    matrix(exp(c(l*log(0.68), l*log(0.68))),
                                           nrow=length(e), ncol=2)
                                },
                                gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                      designs=list(paramMAMS, paramCPBAL, paramCPRAR))

paramMAH1.FS <- list(name="MAH1.FS", backend=survEvd,
                     bpar=list(distr=fspline,
                               taper=function(r, e) {
                                   matrix(exp(c(rep(0, 2*length(e))*(r == 0) +
                                            rep(log(0.42), 2*length(e))*(r != 0))),
                                          nrow=length(e), ncol=2)
                               },
                               gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                     designs=list(paramMAMS, paramCPBAL, paramCPRAR))

paramMAH2.FS <- list(name="MAH2.FS", backend=survEvd,
                     bpar=list(distr=fspline,
                               taper=function(r, e) {
                                   l <- pmin(e/250, 1)                                   
                                   matrix(exp(c(l*log(0.68), l*log(0.68))*(r == 0) +
                                          c((1-l)*log(0.8) + l*log(0.5), (1-l)*log(0.8) + l*log(0.5))*(r == 1) +
                                          c((1-l)*log(0.6) + l*log(0.3), (1-l)*log(0.6) + l*log(0.3))*(r == 2) +
                                          c((1-l)*log(0.42) + l*log(0.2), (1-l)*log(0.42) + l*log(0.2))*(r == 3) +
                                          c((1-l)*log(0.42) + l*log(0.2), (1-l)*log(0.42) + l*log(0.2))*(r == 4)),
                                          nrow=length(e), ncol=2)
                               },
                               gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                     designs=list(paramMAMS, paramCPBAL, paramCPRAR))

#################################################################
#################################################################
#################################################################

param2AH0.NP <- list(backend=survEvd,
                     bpar=list(distr=evd.np, HR.R=c(1, 1), HR.D=c(1, 1),
                               gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                     designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                  paramBayesGSD, paramBGSDRAR))

param2AH1.NP <- list(backend=survEvd,
                     bpar=list(distr=evd.np, HR.R=c(1, 0.46), HR.D=c(1, 0.46),
                               gensurv=genSurv, e.max=10, trigger=10, setup=setup),
                     designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
                                  paramBayesGSD, paramBGSDRAR))

########## SCM

## param2AH0 <- list(backend=scm,
##                   bpar=list(dcm=list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
##                                      beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39,
##                                      delta2=1/9.42, tbegin=0, tend=700, changePoint=224,
##                                      f1=0.52, f2=0.21,
##                                      HR.D1=c(1, 1), HR.D2=c(1, 1)),
##                            func=scmSEIRD, gensurv=genScmSEIRD,
##                            e.max=10, trigger=10, setup=setup), ## antiTrigger=20,
##                   designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
##                                paramBayesGSD, paramBGSDRAR))

## param2AH1 <- list(backend=scm,
##                   bpar=list(dcm=list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
##                                 beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39, delta2=1/9.42,
##                                 tbegin=0, tend=700, changePoint=224, f1=0.52, f2=0.21,
##                                 HR.D1=c(1, 0.5), HR.D2=c(1, 0.5)),
##                            func=scmSEIRD, gensurv=genScmSEIRD,
##                            e.max=10, trigger=10, setup=setup), ## antiTrigger=20,
##                   designs=list(param2ACC, param2AparHC, paramGSD, paramOSBAL,
##                                paramBayesGSD, paramBGSDRAR))

## paramMAH0 <- list(backend=scm,
##                   bpar=list(dcm=list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
##                                 beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39, delta2=1/9.42,
##                                 tbegin=0, tend=700, changePoint=224, f1=0.52, f2=0.21,
##                                 HR.D1=c(1, 1, 1, 1, 1), HR.D2=c(1, 1, 1, 1, 1)),
##                            func=scmSEIRD, gensurv=genScmSEIRD,
##                            e.max=10, trigger=10, setup=setup), ## antiTrigger=20,
##                   designs=list(paramMAparCC, paramCPBAL, paramCPRAR))

## paramMAH1 <- list(backend=scm,
##                   bpar=list(dcm=list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
##                                 beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39, delta2=1/9.42,
##                                 tbegin=0, tend=700, changePoint=224, f1=0.52, f2=0.21,
##                                 HR.D1=c(1, rep(0.5, 4)), HR.D2=c(1, rep(0.5, 4))),
##                            func=scmSEIRD, gensurv=genScmSEIRD,
##                            e.max=10, trigger=10, setup=setup), ## antiTrigger=20,
##                   designs=list(paramCPBAL, paramCPRAR))

## paramMAH2 <- list(backend=scm,
##                   bpar=list(dcm=list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
##                                 beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39, delta2=1/9.42,
##                                 tbegin=0, tend=700, changePoint=224, f1=0.52, f2=0.21,
##                                 HR.D1=c(1, 1, 1, 0.75, 0.5), HR.D2=c(1, 1, 1, 0.75, 0.5)),
##                            func=scmSEIRD, gensurv=genScmSEIRD,
##                            e.max=10, trigger=10, setup=setup), ## antiTrigger=20,
##                   designs=list(paramCPBAL, paramCPRAR))

designNames <- c("2ACC", "2AparHC", "GSD", "OS.BAL", "GSD.Bayes", "GSD.Bayes.RAR",
                 "CP.BAL", "CP.RAR", "MAparCC", "MSA", "MAMS", NA)

trialResults <- c("futility", "success", "l.max", "fail", NA)


simTrial <- function(trial.id, param) {
    x <- param$backend(param$bpar)

    df0 <- summary(x)
    df0$design <- "0"
    df0$trial.id <- trial.id

    if(nrow(x$data1) > 10) {
        ## simulate each trial design in param$designs
        df1 <- do.call(rbind, lapply(param$designs, function(el) {
            el$func(x$data1, x$datah, el$param, param$bpar)
        }))
        df1 <- as.data.frame(df1)
    } else { ## trial program never triggered
        df1 <- data.frame(design=length(designNames), result=length(trialResults),
                          arm=NA, ia=NA, time=NA, n=NA, p=NA, su=NA, fu=NA, dprev=NA)
    }

    df1$design <- designNames[df1$design]
    df1$result <- trialResults[df1$result + 1]
    df1$duration <- df1$time - df0$t1
    df1$trial.id <- trial.id
    
    list(df0=df0, df1=df1)   
}

sim <- function(M, param, src.file="sim.R", trial.fun=simTrial, n.cores=1, seed=NULL) {
    f <- function(i) trial.fun(i, param)
    
    if(n.cores > 1) {
        require(parallel)
        
        cl <- makeCluster(n.cores)

        if(!is.null(seed)) clusterSetRNGStream(cl, seed)

        clusterExport(cl,
                      c("src.file", "param", "evd.data", "fspline", "evd.np"),
                      envir=environment())
        
        clusterEvalQ(cl, {
            source(src.file)
            NULL
        })
        
        on.exit(stopCluster(cl))
        
        x <- parLapply(cl, 1:M, f)
    } else {
        if(!is.null(seed)) set.seed(seed)
        x <- lapply(1:M, f)
    }
    
    df0 <- do.call(rbind, lapply(x, function(xx) xx[[1]]))
    df1 <- do.call(rbind, lapply(x, function(xx) xx[[2]]))
    x <- list(df0=df0, df1=df1)
    
    class(x) <- c("sim", class(x))
    x
}

superSim2 <- function(M, n.cores=1) {
    lst <- list(paramMAH1.FS, paramMAH2.FS)
    lapply(lst, function(el) sim(M, el, n.cores=n.cores))
}

print.sim <- function(x, ...) {
    print(x$df0 %>%
          filter(is.finite(t0) & is.finite(t1)) %>%
          summarise(min.t0=min(t0), mean.t0=round(mean(t0)), max.t0=max(t0),
                    min.t1=min(t1), mean.t1=round(mean(t1)), max.t1=max(t1)))
    print(paste("Trials:", max(x$df1$trial.id)))
    print(paste("Designs:", do.call(paste, as.list(unique(x$df1$design)))))
}

summary.sim <- function(x, ...) {
    df0 <- x$df0 %>%
        filter(is.finite(t0) & is.finite(t1)) %>%
        summarise(min.t0=min(t0), mean.t0=round(mean(t0)), max.t0=max(t0),
                  min.t1=min(t1), mean.t1=round(mean(t1)), max.t1=max(t1),
                  min.n0=min(n0), mean.n0=mean(n0), max.n0=max(n0),
                  min.n1=min(n1), mean.n1=mean(n1), max.n1=max(n1))
    
    df0$pt0 <- mean(!is.finite(x$df0$t0))
    df0$pt1 <- mean(!is.finite(x$df0$t1))

    df1 <- x$df1 %>%
        group_by(design, trial.id) %>%
        filter(!is.na(design)) %>%
        summarise(n=max(n), time=max(time), duration=max(duration),
                  su=any(su == 1, na.rm=TRUE), suNone=all(su == 0, na.rm=TRUE),
                  fu=any(fu == 1, na.rm=TRUE), p=mean(p)) %>%
        summarise(suAny=mean(su), suNone=mean(suNone), fuAny=mean(fu), p=mean(p),
                  mean.n=mean(n), max.n=max(n),
                  mean.time=mean(time), max.time=max(time), min.time=min(time),
                  mean.duration=mean(duration), max.duration=max(duration))

    df2 <- x$df1 %>%
        group_by(design, arm, trial.id) %>%
        filter(!is.na(design)) %>%
        summarise(n=max(n), time=max(time), duration=max(duration), su=any(su == 1, na.rm=TRUE),
                  fu=any(fu == 1, na.rm=TRUE), p=mean(p)) %>%
        summarise(su=mean(su), fu=mean(fu), p=mean(p),
                  min.n=min(n), mean.n=mean(n), max.n=max(n),
                  time=mean(time), mean.duration=mean(duration), max.duration=max(duration))

    df3 <- x$df1 %>%
        group_by(design, ia, arm, trial.id) %>%
        filter(!is.na(design)) %>%
        summarise(n=max(n), time=max(time), duration=max(duration), su=any(su == 1, na.rm=TRUE),
                  fu=any(fu == 1, na.rm=TRUE), p=mean(p)) %>%
        summarise(su=mean(su), fu=mean(fu), p=mean(p), min.n=min(n), mean.n=mean(n), max.n=max(n),
                  time=mean(time), mean.duration=mean(duration), max.duration=max(duration))

    df4 <- x$df1 %>% group_by(design, result, ia, arm) %>%
        summarise(n=n()) %>%
        group_by(design, result) %>%
        summarise(n=max(n))

    df5 <- x$df1 %>%
        group_by(design, trial.id, ia) %>%
        summarise(n=sum(n)) %>%
        summarise(n=max(n)) %>%
        summarise(n=mean(n))

    df6 <- x$df1 %>%
        group_by(design, trial.id) %>% summarise(dprev=max(dprev)) %>%
        summarise(dprev=mean(dprev))

    if(length(unique(x$df1$arm)) > 2) {
        df7 <- x$df1 %>% group_by(design, trial.id, arm) %>% filter(!is.na(design)) %>% summarise(su=any(su==1,na.rm=TRUE)) %>% summarise(su=any(su[3:4])) %>% summarise(su=mean(su))
    } else df7 <- NULL
    
    list(df0=df0, df1=df1, df2=df2, df3=df3, df4=df4, df5=df5, df6=df6, df7=df7)
}

plotTimeCourse <- function(tc=get(load("timecourse.res"))) {
    ##tc <- getTimeCourse(data)
    ##plot(tc$time, tc$total, type="l", xlab="Time [days]", ylab="Total number conf. + hosp.")
    p <- ggplot(tc, aes(x=time, y=total)) + geom_line() + xlab("Time [days]") + ylab("Total number confirmed and hospitalized") + theme_bw()
    
    ##abline(v=df$df0$mean.t0, col="green")
    ##p <- p + geom_vline(aes(xintercept=11), linetype=2)
    ##abline(v=df$df0$min.t0, lty=2)
    ##abline(v=df$df0$max.t0, lty=2)
    
    ##abline(v=df$df0$mean.t1, col="blue")
    p <- p + geom_vline(aes(xintercept=53), linetype=2)
    p <- p + geom_vline(aes(xintercept=111), linetype=2)
    p <- p + geom_vline(aes(xintercept=211), linetype=2)
    p <- p + geom_vline(aes(xintercept=311), linetype=2)
    p 
}
