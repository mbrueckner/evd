#' interim.gsd
#'
#' implementation of param$interim interface
#' 
#' @param data data.frame of interim data
#' @param param list of design specific parameters
#' @param l number of current interim analysis
#' @param code result code of previous stage (not used in GSD design)
#' @return data.frame
interim.gsd <- function(data, param, l, code) {
    ed <- getEndpoint(data, param$x, TRUE)

    if(ed$var[1] + ed$var[2] == 0) Z <- 0
    else {
        ##lo <- logodds(ed$n, ed$p, ed$var)
        ##if(!is.finite(lo[1])) Z <- lo[1]        
        ##else Z <- lo[1] / sqrt(lo[2])
        Z <- (ed$p[1] - ed$p[2]) / sqrt(ed$var[1] + ed$var[2])
    }
    
    su <- c(FALSE, Z > param$upper[l])
    fu <- c(FALSE, Z < param$lower[l])

    data.frame(arm=ed$R, n=ed$n, y=ed$y, su=su, fu=fu, p=ed$p,
               new.rp=rep(0.5, 2), theta=c(0, 0),
               n.cur=rep.int(nrow(data), 2), code=rep.int(0, 2),
               ia.time=rep(max(data$entry + data$time), 2))
}

delta <- function(p0, p1) (p0-p1) / sqrt(p0*(1-p0) + p1*(1-p1))

## sample size for fixed 
size <- function(p0, p1, alpha=0.025, beta=0.1) 2 * (qnorm(alpha) + qnorm(beta))^2 / delta(p0, p1)^2

paramGSD <- list(func=driver, param=list(trt.ids=0:1, interim=interim.gsd, design=3,
                                         alpha=0.025, anaDelay=1, l.max=5, x=14,
                                         n.max=220, c.max=47, reallocate=FALSE))

require(gsDesign)
## non-binding futility stopping
tmp <- gsDesign(k=paramGSD$param$l.max, test.type=3, beta=0.1, delta=delta(0.40, 0.20))

paramGSD$param$n.max <- ceiling(max(tmp$n.I*2))
paramGSD$param$c.max <- floor(paramGSD$param$n.max / paramGSD$param$l.max)
paramGSD$param$upper <- tmp$upper$bound
paramGSD$param$lower <- tmp$lower$bound ##c(rep(0, paramGSD$param$l.max-1), paramGSD$param$upper[paramGSD$param$l.max])

## no futility stopping for H0 sim
paramGSD0 <- paramGSD
paramGSD0$param$lower <- rep(-Inf, length(paramGSD0$param$lower))

param2ACC <- paramGSD
param2ACC$param$design <- 1
param2ACC$param$l.max <- 1
param2ACC$param$n.max <- ceiling(size(0.4, 0.2))
param2ACC$param$c.max <- ceiling(size(0.4, 0.2))
param2ACC$param$ifrac <- 1
param2ACC$param$lower <- qnorm(0.975)
param2ACC$param$upper <- qnorm(0.975)
    
