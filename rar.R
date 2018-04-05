## Closed Platform Response-Adaptive Randomization Design

pmf.beta.binomial <- function(n, k, a, b) choose(n, k) * beta(k + a, n - k + b) / beta(a, b)

postParam <- function(n, y, a, b) c(a + y, b + n - y)

postProb <- function(delta, n, y, priorA, priorB) {
    ## calculate parameters of posterior beta distribution
    postA <- c(priorA[1] + y[1], priorA[2] + n[1] - y[1])
    postB <- c(priorB[1] + y[2], priorB[2] + n[2] - y[2])

    ## P(p_A > p_B + delta | nA, nB, yA, yB)
    f <- function(x) pbeta(x - delta, postA[1], postA[2]) * dbeta(x, postB[1], postB[2])
    integrate(f, 0, 1)$value
}

## success criterion: P(p_i < p_0 | data) > S
checkSuccess <- function(S, n, y, prior) {
    x <- sapply(2:length(n), function(i) postProb(0, n[c(i, 1)], y[c(i, 1)], prior[i,], prior[,1]))
    c(FALSE, x > S)
}

## futility criterion: P(p_i < p_0 - delta | data) < F
checkFutility <- function(F, delta, n, y, prior, n.max) {
    x <- sapply(2:length(n), function(i) postProb(delta, n[c(i, 1)], y[c(i, 1)], prior[i,], prior[1,]))
    c(FALSE, x < F)
}

popt <- function(alpha, beta) {
    f <- function(x, i) prod(1 - pbeta(x, alpha[-i], beta[-i]))*dbeta(x, alpha[i], beta[i])
    vf <- Vectorize(f, vectorize.args="x")
    sapply(1:length(alpha), function(i) integrate(vf, lower=0, upper=1, i=i)$value)
}

## balanced allocation rule
adjustRP.dummy <- function(alpha, beta, theta, n, fu) {
    rp <- rep.int(0, length(n))
    m <- sum(fu == 0)
    rp[fu == 0] <- rep(1/m, m)
    rp
}

## Response-adaptive randomization as in Saville/Berry (2016)
adjustRP.SB <- function(alpha, beta, theta, n, fu) {
    rp <- rep.int(0, length(n))
    varP <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))                
    rp[fu == 0] <- sqrt(theta[fu == 0] * varP[fu == 0] / (n[fu == 0] + 1))
    rp / sum(rp)
}

#' interim.bayes
#'
#' Implementation of param$interim interface
#' 
#' @param data data.frame of interim data
#' @param param list of design specific parameters
#' @param l number of current interim analysis
#' @param code result code of previous stage (not used)
#' @return data.frame
interim.bayes <- function(data, param, l, code) {
    delta <- param$delta
    S <- param$S
    F <- param$F
    n.max <- param$n.max

    ed <- getEndpoint(data, param$x, FALSE)
    n <- ed$n
    y <- ed$y

    trt.ids <- unique(data$R)
    prior <- param$prior[trt.ids+1,]
    
    su <- checkSuccess(param$S, n, y, prior)
    fu <- checkFutility(param$F, param$delta, n, y, prior, param$n.max)
    fu[su] <- 0 # happens when success criterion TRUE and n.max reached at the same time

    ## posterior parameters
    alpha <- prior[,1] + y
    beta <- prior[,2] + n - y
    theta <- popt(alpha, beta)

    ## posterior mode
    mode <- (alpha - 1) / (alpha + beta - 2)
    
    ## recalculate randomization probabilities
    if(any(fu == 0)) {
        rp <- param$adjustRP(alpha, beta, theta, n, fu)
    }

    data.frame(arm=ed$R, n=n, y=y, su=su, fu=fu,
               p=mode, new.rp=rp, theta=theta,
               n.cur=rep.int(nrow(data), length(n)),
               code=rep.int(0, length(n)), ia.time=rep(max(data$entry + data$time), length(n)))
}

## closed-platform response-adaptive randomization (multi-arm, multi-stage)
paramCPRAR <- list(name="CP.RAR", func=driver,
                   param=list(l.max=5, design=8, anaDelay=1,
                              interim=interim.bayes,
                              adjustRP=adjustRP.SB,
                              trt.ids=0:4,
                              prior=matrix(1, nrow=5, ncol=2),
                              c.max=110, ## cohort size
                              n.max=550, ## group maximum
                              delta=0.1,
                              S=0.99, ## success boundary
                              F=0.1, x=14, ## futility boundary
                              reallocate=FALSE))

## no futility stopping for H0 sim
paramCPRAR0 <- paramCPRAR
paramCPRAR0$param$F <- 0

## closed-platform balanced allocation (multi-arm multi-stage)
paramCPBAL <- list(name="CP.BAL", func=driver,
                      param=list(l.max=5, design=7, anaDelay=1,
                                 interim=interim.bayes,
                                 adjustRP=adjustRP.dummy,
                                 trt.ids=0:4,
                                 prior=matrix(1, nrow=5, ncol=2),
                                 c.max=110, ## cohort size
                                 n.max=550, ## group maximum
                                 delta=0.1,
                                 S=0.99, ## success boundary
                                 F=0.1, x=14, ## futility boundary
                                 reallocate=FALSE))

## no futility stopping for H0 sim
paramCPBAL0 <- paramCPBAL
paramCPBAL0$param$F <- 0

## closed-platform response-adaptive randomization
paramBGSDRAR <- list(name="GSD.Bayes.RAR", func=driver,
                   param=list(l.max=5, K=1, design=6, anaDelay=1,               
                              interim=interim.bayes,
                              adjustRP=adjustRP.SB,
                              trt.ids=0:1,
                              prior=matrix(1, nrow=2, ncol=2),
                              c.max=43, ## cohort size
                              n.max=216, ## group maximum
                              delta=0.1,
                              S=0.99, ## success boundary
                              F=0.1, x=14, ## futility boundary
                              reallocate=FALSE))

## no futility stopping for H0 sim
paramBGSDRAR0 <- paramBGSDRAR
paramBGSDRAR0$param$F <- 0

## two-arm Bayesian GSD
paramBayesGSD <- list(name="GSD.Bayes", func=driver,
                      param=list(l.max=5, K=1, design=5, anaDelay=1,
                                 interim=interim.bayes,
                                 adjustRP=adjustRP.dummy,
                                 trt.ids=0:1,
                                 prior=matrix(1, nrow=2, ncol=2),
                                 c.max=43, ## cohort size
                                 n.max=216, ## group maximum
                                 delta=0.1,
                                 S=0.99, ## success boundary
                                 F=0.1, x=14, ## futility boundary
                                 reallocate=FALSE))

## no futility stopping for H0 sim
paramBayesGSD0 <- paramBayesGSD
paramBayesGSD0$param$F <- 0

## single-stage balanced allocation
paramOSBAL <- list(name="OS.BAL", func=driver,
                   param=list(l.max=1, K=1, design=4, anaDelay=1,
                              interim=interim.bayes,
                              adjustRP=adjustRP.dummy,
                              trt.ids=0:1,
                              prior=matrix(1, nrow=2, ncol=2),
                              c.max=211, ## cohort size
                              n.max=211, ## group maximum
                              delta=0.1,
                              S=0.975, ## success boundary
                              F=0.1, x=14, ## futility boundary
                              reallocate=FALSE))

## no futility stopping for H0 sim
paramOSBAL0 <- paramOSBAL
paramOSBAL0$param$F <- 0
