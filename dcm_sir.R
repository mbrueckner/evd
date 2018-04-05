source("import.R")
source("data.R")

library(deSolve)


gradSEIR <- function(t, x, vpar) {
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]

    gamma <- vpar$gamma    
    beta <- vpar$beta
    delta <- vpar$delta
    
    N <- S + E + I + R
    dS <- -beta(t)*S*I/N
    dE <- beta(t)*S*I/N - gamma*E
    dI <- gamma*E - delta(t)*I
    dR <- delta(t)*I
    
    list(c(dS, dE, dI, dR))
}

gradSEIRD <- function(t, x, vpar) {
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]
    D <- x[5]

    gamma <- vpar$gamma
    
    if(t < vpar$changePoint) {
        beta <- vpar$beta1
        delta <- vpar$delta1
        f <- vpar$f1
    } else {
        beta <- vpar$beta2
        delta <- vpar$delta2
        f <- vpar$f2
    }
    
    N <- S + E + I + D + R
    dS <- -beta*S*I/N
    dE <- beta*S*I/N - gamma*E
    dI <- gamma*E - delta*I
    dR <- (1-f)*delta*I
    dD <- f*delta*I
    
    list(c(dS, dE, dI, dR, dD))
}

paramSL <- list(N=6.316e6, I0=1, E0=8, R0=0, D0=0,
                beta1=0.15, beta2=0.087, gamma=1/10, delta1=1/10.39, delta2=1/9.42,
                tbegin=0, tend=700, changePoint=224, f1=0.52, f2=0.21)

paramEVD <- list(N=23.5e6, I0=100, E0=5000, R0=0, D0=0,
                 gamma=1/10,
                 tbegin=0, tend=700,
                 beta=na$h1, delta=na$h2)

## beta = transmission rate [1/Nt]
## gamma = infectious rate [1/t], 1/gamma = average incubation period
## delta = recovery/death rate [1/t], 1/delta = average infectious period = average time to recovery/death

SEIRD <- function(param) {    
    df <- as.data.frame(with(param, {
        vt <- seq(tbegin, tend, 1)  
        vpar <- list(beta1=beta1, beta2=beta2, gamma=gamma, delta1=delta1, delta2=delta2,
                     changePoint=changePoint, f1=f1, f2=f2)
        inits <- c(S=N-I0-E0-R0, E=E0, I=I0, R=0, D=0)  ## order must match model!
        lsoda(inits, vt, gradSEIRD, vpar)
    }))

    ## number of newly infected patients
    df$newI <- c(1, round(-diff(df$S) / diff(df$time), 0))
    
    df <- round(df, 0)    
    df
}

SEIR <- function(param) {    
    df <- as.data.frame(with(param, {
        vt <- seq(tbegin, tend, 1)  
        vpar <- list(gamma=gamma, beta=beta, delta=delta)
        inits <- c(S=N-I0-E0-R0, E=E0, I=I0, R=0)  ## order must match model!
        lsoda(inits, vt, gradSEIR, vpar)
    }))

    ## number of newly infected patients
    df$newI <- c(1, round(-diff(df$S) / diff(df$time), 0))    
    df <- round(df, 0)    
    df
}

plotDCM <- function(param) {
    sl <- ebolaSL() ## Sierra Leone data
    x <- SEIRD(param)

    df <- data.frame(time=c(sl$days, x$time), type=c(rep("Data", nrow(sl)), rep("Model", nrow(x))), inc=c(sl$inc, cumsum(x$newI)), deaths=c(sl$death, x$D))

    ggplot(df, aes(x=time, y=inc, linetype=type, col="Infections")) + geom_line(size=1) +
        geom_line(aes(x=time, y=deaths, linetype=type, col="Deaths"), size=1) + theme_bw() +
        scale_color_manual(name = "Variable", values = c(Infections = I("black"), Deaths = I("red"))) +
        scale_linetype_manual(name = "Source", values = c(Data="solid", Model="dotted")) +
        labs(x="Time", y="Count") + theme(axis.text=element_text(size=13),
                                          axis.title=element_text(size=14),
                                          legend.text=element_text(size=13),
                                          legend.title=element_text(size=14))                                             
}

SEIRD.opt <- function() {
    sl <- ebolaSL()

    require(nloptr)
    
    f <- function(x) {
        param$beta1 <- x[1]
        param$beta2 <- x[2]
        param$delta1 <- 1/x[3]
        param$delta2 <- 1/x[4]
        param$changePoint <- x[5]
        param$f1 <- x[6]
        param$f2 <- x[7]
        param$E0 <- x[8]
        
        df <- SEIRD(param)

        gi <- approxfun(df$time, cumsum(df$newI), method="linear")
        gd <- approxfun(df$time, df$D, method="linear")
        
        log(sum((gi(sl$days) - sl$inc)^2 + (gd(sl$days) - sl$death)^2))
    }

    fit0 <- nloptr(c(0.15, 0.087, 10.39, 9.41, 224, 0.52, 0.21, 8), f, ##c(0.15, 0.087, 10.38, 9.37, 225, 0.53, 0.2, 8), f,
                   opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=0.00001, maxeval=150))
}
