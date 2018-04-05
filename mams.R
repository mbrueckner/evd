#' interim.mams
#'
#' Implementation of param$interim interface for MAMS trial
#' 
#' @param data data.frame of interim data
#' @param param list of design specific parameters
#' @param l number of current interim analysis
#' @param code result code of previous stage (not used in GSD design)
#' @return data.frame
interim.mams <- function(data, param, l, code) {
    ed <- getEndpoint(data, param$x, TRUE)
    
    D <- nrow(ed) - 1

    if(ed$var[1] == 0) Z <- rep(0, D)
    else {
        tmp <- sapply(1:D, function(i) {
            ##lo <- logodds(ed$n[c(1,i)], ed$p[c(1,i)], ed$var[c(1,i)])
            ##if(!is.finite(lo[1])) lo[1]
            ##else lo[1] / lo[2]
           
            if(ed$var[i+1] == 0) 0
            else ed$p[1] - ed$p[i+1]
        })
        
        Im <- ed$n[1] / (0.4*0.6*2)
        Imax <- param$Imax
        Z <- tmp * Im / sqrt(Imax)
    }
    
    su <- c(FALSE, Z > param$upper[l])
    fu <- c(FALSE, Z < param$lower[l])

    new.rp <- rep(1/sum(!fu), D+1)
    new.rp[fu] <- 0

    data.frame(arm=ed$R, n=ed$n, y=ed$y, su=su, fu=fu, p=ed$p,
               new.rp=new.rp, theta=rep(0, D+1),
               n.cur=rep.int(nrow(data), D+1), code=rep.int(0, D+1),
               ia.time=rep(max(data$entry + data$time), D+1))
}

paramMAMS <- list(func=driver, param=list(trt.ids=0:4, interim=interim.mams, design=11,
                                          ifrac=seq(0.2, 1, by=0.2), alpha=0.025, anaDelay=1,
                                          l.max=5, x=14, n.max=600, c.max=120,
                                          reallocate=FALSE))

##tmp <- gsDesign(k=paramMAMS$param$l.max, test.type=3, beta=0.1, n.fix=paramGSD$param$n.max/2)

D <- length(paramMAMS$param$trt.ids)-1
J <- paramMAMS$param$l.max

##b <- bnd(D=D, J=J, alpha=cumsum(tmp$upper$spend), ##beta=cumsum(tmp$lower$spend),
##         sigma=rep(1, D+1), alloc=rep(1, D+1), N=1, n0=tmp$n.I, lattice=mcLattice) 
              
paramMAMS$param$n.max <- 550 ##ceiling(max(tmp$n.I*2))
paramMAMS$param$c.max <- 110 ##ceiling(paramMAMS$param$n.max / paramMAMS$param$l.max)
paramMAMS$param$Imax <- 110 / (0.4*0.6*2)
paramMAMS$param$upper <- c(3.628937, 3.391564, 3.039619, 2.806660, 2.539892)
paramMAMS$param$lower <- c(-0.67072001, -0.40570110, -0.04258631,  0.33013320, paramMAMS$param$upper[5])

paramMAMS0 <- paramMAMS
paramMAMS0$param$lower <- c(rep(-Inf, 4), paramMAMS$param$lower[5])
                                
