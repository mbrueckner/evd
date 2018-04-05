
getRawData <- function() {
    data <- read.csv("evd.csv", sep=",")
    
    data$DateOnsetInferred <- as.Date(data$DateOnsetInferred, "%d/%m/%Y")
    data$DateOutcomeComp <- as.Date(data$DateOutcomeComp, "%d/%m/%Y")
    data$DateHospitalCurrentAdmit <- as.Date(data$DateHospitalCurrentAdmit, "%d/%m/%Y")

    data <- data[!is.na(data$DateOnsetInferred), ]
    data
}

getSIRData <- function(file="evd.csv") {
    data0 <- read.csv(file, sep=",")
    data <- data0[!is.na(data0$DateOnsetInferred) & !is.na(data0$DateOutcomeComp), ]
    
    data$DateOnsetInferred <- as.Date(data$DateOnsetInferred, "%d/%m/%Y")
    data$DateOutcomeComp <- as.Date(data$DateOutcomeComp, "%d/%m/%Y")    
    
    origin <- min(data$DateOnsetInferred)
    data$I <- as.numeric(difftime(data$DateOnsetInferred, origin, units="days"))
    data$R <- as.numeric(difftime(data$DateOutcomeComp, origin, units="days"))
    
    data <- data[order(data$I), ]
    data$id <- 1:nrow(data)    
    data[data$I <= data$R, ]
}

getNA <- function(data) {
    data$from <- 0
    data$to <- 1
    data$time <- data$I
    data$start <- 0
    data$stop <- data$time
    data$status <- 1
    fit1 <- mvna(data, state.names=c("0", "1"), tra=tra, cens.name="cens")
    fs1 <- flexsurvspline(Surv(time, status) ~ 1, data=data, k=6)
    
    data$from <- 0
    data$to <- 2
    data$time <- data$R
    data$start <- data$I
    data$stop <- data$R
    
    fit2 <- mvna(data, state.names=c("0", "2"), tra=tra, cens.name="cens")
    fs2 <- flexsurvspline(Surv(time, status) ~ 1, data=data, k=6)

    H1 <- Hsurvspline(fit1[[1]]$time, fs1$opt$par, knots=fs1$knots)
    H2 <- Hsurvspline(fit2[[1]]$time, fs2$opt$par, knots=fs2$knots)

    f1 <- function(u) hsurvspline(u, fs1$opt$par, knots=fs1$knots)
    f2 <- function(u) hsurvspline(u, fs2$opt$par, knots=fs2$knots)
    
    list(h1=f1, h2=f2)
    
    ##data.frame(time=c(fit1[[1]]$time, fit2[[1]]$time),
    ##           na=c(fit1[[1]]$na, fit2[[1]]$na, H1, H2),
    ##           state=c(rep("I", length(fit1[[1]]$time)), rep("R", length(fit2[[1]]$time))),
    ##           method=rep(c("na", "spline"), each=length(fit1[[1]]$time)+length(fit2[[1]]$time)))
               
}

getEVDData <- function(output=FALSE) {
    ## Patient-level EVD data
    data0 <- read.csv("evd.csv", sep=",")
    
    ## confirmed and hospitalized cases only
    data <- data0[(data0$HospitalizedEver == "Yes") & (data0$EpiCaseDef == "confirmed"), ]

    ## remove missing data
    data <- data[!is.na(data$DateOutcomeComp) &
                 !is.na(data$DateHospitalCurrentAdmit) &
                 !is.na(data$FinalStatus), ]

    data$DateOnsetInferred <- as.Date(data$DateOnsetInferred, "%d/%m/%Y")
    data$DateOutcomeComp <- as.Date(data$DateOutcomeComp, "%d/%m/%Y")
    data$DateHospitalCurrentAdmit <- as.Date(data$DateHospitalCurrentAdmit, "%d/%m/%Y")

    ## hospitalization times == (potential) entry times
    origin <- min(data$DateHospitalCurrentAdmit)
    data$entry <- as.numeric(difftime(data$DateHospitalCurrentAdmit, origin, units="days"))

    ## days since infection
    data$time0 <- as.numeric(difftime(data$DateOutcomeComp, data$DateOnsetInferred, units="days"))

    ## days since entry
    data$time1 <- as.numeric(difftime(data$DateOutcomeComp, data$DateHospitalCurrentAdmit, units="days"))

    ## patients being infected and hospitalized on the same day
    data$time0[data$time0 == 0] <- 0.25

    ## patients being hospitalized and dying/recovering on the same day
    data$time1[data$time1 == 0] <- 0.5

    data$status <- rep("cens", nrow(data))
    data$status[data$FinalStatus == "Alive"] <- 1  ## recovery
    data$status[data$FinalStatus == "Dead"] <- 2  ## death

    data$from <- 0
    data$to <- data$status
    
    data$d14 <- (data$status == 2) & (data$time1 <= 14)
    data$d14 <- (data$status == 2) & (data$time1 <= 28)

    data$time <- data$time1
    
    data$centre <- data$TreatmentCentre
    data$country <- data$Country

    data$id <- 1:nrow(data)
    
    evd.data <- data[order(data$entry), ]

    if(output) save(evd.data, file="evd.RData")          

    evd.data
}

getTimeCourse <- function(data) {
    origin2 <- min(data$DateOnsetInferred)

    df <- data.frame(time0=as.numeric(difftime(data$DateOnsetInferred, origin2, units="days")), d14=data$d14)
                     ##time1=as.numeric(difftime(data$DateOutcomeComp, origin2, units="days")))

    df1 <- df %>% group_by(time0) %>% summarise(n=n(), new=n(), death=sum(d14)) %>% arrange(time0)
    ##df2 <- df %>% group_by(time1) %>% summarise(n=-n()) %>% arrange(time1)
    
    colnames(df1) <- c("time", "n", "new", "death")
    ##colnames(df2) <- c("time", "n")
    
    ##df <- full_join(df1, df2) %>% group_by(time) %>% summarise(n=sum(n), new=max(new))
    ##df$new[is.na(df$new)] <- 0
    df1$n <- cumsum(df1$n)
    df1$total <- cumsum(df1$new)
    df1$death <- cumsum(df1$death)
    df1
}

getCMA <- function(data, start=0) {
    df <- data.frame(time=data$entry, p=cumsum(data$d14) / 1:nrow(data))
    df[df$time >= start,]
}

getRates <- function(data, t0=0) {
    data.frame(time=t0,
               n=sapply(t0, function(start) sum(data$entry >= start)),
               p=sapply(t0, function(start) mean(data$d14[data$entry >= start])))
}
