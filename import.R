ebola <- function() {
    require(dplyr)
    
    data <- getSIRData("evd.csv")
    
    df1 <- data %>% group_by(I) %>% summarise(n=n())
    df2 <- data %>% group_by(R) %>% summarise(n=n())

    colnames(df1) <- c("days", "I")
    colnames(df2) <- c("days", "R")

    df <- full_join(df1, df2)

    df$I[is.na(df$I)] <- 0
    df$R[is.na(df$R)] <- 0

    df <- df[order(df$days),]

    df$cI <- cumsum(df$I)
    df$cR <- cumsum(df$R)

    df
}
