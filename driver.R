#' Determine calendar time at which a certain sample size / event number has been reached
#'
#' @param data data.frame
#' @param n number of patients (or number of events if events=TRUE)
#' @param x number of days
#' @param events boolean variable
#' @return vector of length two (first element=interim analysis time, second element: TRUE if requested number of events/patients reached and FALSE otherwise)
#' 
getInterimTime <- function(data, x=14, target=nrow(data), events=FALSE) {
    n <- nrow(data)
    Y <- data$entry + pmin(data$time, x)
    Y <- sort(Y)
    
    if(events) {
        stop("TODO")
        Y <- Y[data$status == 2 & data$time <= x]
        ia.time <- Y[min(length(Y), n)]
    } else {
        ia.time <- Y[min(length(Y), n)]
    }
    c(ia.time, length(Y) < target)
}

## transition matrix for competing risks model to be used by \link{\code{etm}}
cif.tra <- matrix(FALSE, nrow=3, ncol=3)
cif.tra[1,2:3] <- TRUE

#' Estimate and evaluate CIF at Day x
#' 
#' @param data data.frame containing variables id, time, from, to as required by \link{\code{etm}}
#' @param x time [days] to evaluate the CIF at
#' @return 2-vector of CIF value at x and the corresponding variance
#' 
cifx <- function(data, x) {
    res <- etm(data, state.names=c("0", "1", "2"), tra=cif.tra, cens.name="cens",
               s=0, t=x, covariance=TRUE)
    
    v <- res$cov[7, 7, dim(res$cov)[3]]
    
    if(is.null(v)) v <- 0
        
    c(res$est[1,3,dim(res$est)[3]], v)
}

#' getEndpoint
#'
#' @param data data.frame
#' @param x number of days
#' @param est.cif if TRUE estimate death rate and its variance by evaluating the CIF
#' @return data.frame with number of patients (n), hits (y), estimated death rate (p) and estimated variance of p (var) per arm
#' 
getEndpoint <- function(data, x, est.cif=TRUE) {
    tmp <- do.call(rbind, lapply(sort(unique(data$R)), function(r) {
        dd <- data[data$R == r, ]
        n <- nrow(dd)
        y <- sum(dd$status == 2 & dd$time <= x)
        st <- c(r, n, y)
        if(est.cif) cf <- cifx(dd, x)
        else {
            p <- y/n
            cf <- c(p, p*(1-p)/n)
        }
        
        c(st, cf)
    }))

    df <- as.data.frame(tmp)
    colnames(df) <- c("R", "n", "y", "p", "var")
    if(any(is.na(df))) browser()
    df
}

#' Calculate log-odds ratio estimator and variance
#'
#' @param n pair of sample sizes
#' @param p pair of estimated probabilities
#' @param v pair of estimated variances of \code{p}
#' @return vector with estimated log-odds ratio and its estimated variance
logodds <- function(n, p, v) c(log(p[1]*(1-p[2])/(p[2]*(1-p[1]))), sum(v / (p^2 * (1-p)^2)))

#' selectData
#'
#' time1 and time2 are in calendar time!
#' patient-wise splitting:
#'       stage1: time1=0, time2=final.time, entry.start=0, entry.stop=ia.time
#'       stage2: time1=0, time2=final.time, entry.start=ia.time, entry.stop=final.time
#' stage-wise splitting:
#'       stage1: time1=0, time2=ia.time, entry.start=0, entry.stop=ia.time
#'       stage2: time1=ia.time, time2=final.time, entry.start=0, entry.limit=final.time
#'
#' @param data data.frame as returned by \link{generateData}
#' @param time1 calendar time (e.g. time of previous interim analysis or 0)
#' @param time2 calendar time (e.g. time of current interim analysis / final time)
#' @param entry.start calendar time, exclude patients recruited before entry.start
#' @param entry.stop calendar time, exclude patients recruited after entry.stop
#'
#' @return data.frame containing only patients recruited before and right-censored at calendar time t
#' 
selectData <- function(data, time1=0, time2=Inf, entry.start=0, entry.stop=Inf) {
    entry.stop <- min(time2, entry.stop)
    data <- data[(data$entry <= entry.stop) & (data$entry >= entry.start) & (data$entry + data$time > time1), ]

    ## maximum of max.n patients
    ##data <- data[1:min(nrow(data), max.n),]

    ## observations with (entry time < start time) are left-truncated (if time1=0 then (time1 - R < 0) thus V=0)
    ##data$V <- pmax(time1 - data$R, 0)

    ## observations with (entry time + event time > time2) are right-censored
    data$status[data$time + data$entry > time2] <- "cens"
    data$to <- data$status
    data$id <- 1:nrow(data)
    
    data$time <- pmin(data$time, pmax(time2 - data$entry, 0.5))
    
    data
}

#' driver
#'
#' Generic sequential trial simulation
#' 
#' @param data data.frame containing all potential entry times
#' @param datah data.frame containing historical control data
#' @param param list of design parameters
#' @param bpar list of backend parameters
#' @return matrix containing "design", "arm", "ia", "time", "n", "p", "su", "fu" variables for each treatment group
driver <- function(data, datah, param, bpar) {
    n.max <- min(param$n.max, nrow(data)) ## maximum group size
    l.max <- param$l.max ## maximum number of interim analyses
    anaDelay <- param$anaDelay ## analysis delay
    
    c.max <- param$c.max ## maximum cohort size
    if(length(c.max) == 1) c.max <- rep(c.max, l.max)
    
    trt.ids <- param$trt.ids
    J0 <- length(param$trt.ids)
    
    ##if(param$historic) data1 <- datah
    data1 <- data.frame() ## accumulating data
    
    res <- NULL
    
    n.cur <- 0
    n.cur0 <- 0
    trialResult <- 0

    ## balanced allocation at beginning
    rp <- rep(1/length(trt.ids), length(trt.ids))

    for(i in 1:l.max) {
        n <- n.max - n.cur

        if(param$reallocate) cm <- c.max[i]
        else cm <- c.max[i] * length(trt.ids) / J0
        
        ##print(c(J0, length(trt.ids), c.max[i], cm))
        
        if(n > 0) { ## get next cohort
            c.cur <- min(n, cm) ## cohort size
            fl <- floor(c.cur*rp)
            d <- c.cur - sum(fl)
           
            if(d > length(trt.ids)) warning("rounding error too large!")

            ## treatment assignment
            R <- c(sample(rep(trt.ids, fl), sum(fl)), sample(trt.ids, d, FALSE))
            
            ## generate survival data for current cohort given treatment assignment
            sdata <- bpar$gensurv(R, data[(n.cur+1):(n.cur+c.cur),], bpar)

            data1 <- rbind(data1, sdata)
            n.cur <- n.cur + nrow(sdata) ## nrow(sdata) == c.cur
        } else { ## target sample size not reached
            trialResult <- 3
            break
        }

        ia.time <- getInterimTime(data1, x=param$x, target=n.cur0+cm, events=FALSE)
        ia.data <- selectData(data1, time1=0, time2=ia.time[1])

        tmp <- param$interim(ia.data, param, i, tmp$code)
        tmp$ia <- i
        tmp$rp <- rp
        tmp$time <- tmp$ia.time + anaDelay
                
        res <- rbind(res, as.matrix(tmp))

        if(is.na(any(tmp$su))) browser()

        ## recruitment target missed, no more patients to continue
        if(ia.time[2] != 0) {
            trialResult <- 3
            break
        }

        ## successful treatment identified
        if(any(tmp$su)) {
            trialResult <- 1
            break
        }
        
        ## all treatments stopped for futility
        if(all(tmp[-1,]$fu != 0)) {
            trialResult <- 0
            break 
        }

        ## last analysis and no success -> futility
        if(i == l.max) {
            trialResult <- 2 ## l.max reached
            break
        }

        trt.ids <- trt.ids[tmp$fu == 0]
        rp <- tmp$new.rp[tmp$fu == 0]
        
        ## remove dropped treatments from data set
        data1 <- data1[data1$R %in% trt.ids,]
        n.cur0 <- nrow(data1)
    }

    if(trialResult == 1) {
        index1 <- which(tmp$su)
        index2 <- which.min(tmp$p[index1])
        best.trt <- trt.ids[index1][index2]

        dprev <- 0 ##deaths.prevented(param$x, max(tmp$ia.time),
                   ##               bpar$HR.D[best.trt+1], bpar$HR.R[best.trt+1],
                   ##               cumhaz=bpar$distr$cumhaz)
    } else dprev <- 0

    res <- cbind(res, dprev=rep(dprev, nrow(res)))
    res <- res[, c("arm", "ia", "time", "n", "p", "su", "fu", "dprev")]
    
    cbind(design=param$design, result=trialResult, res)
}
