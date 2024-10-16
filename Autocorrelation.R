##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: time series qPCR and dPCR measures of sense & antisense RNA.
# Object: autocorrelation times calculated based on qPCR and dPCR measures of sense & antisense RNA. Final output file is used to construct a noise space.
##############################################################

# R functions

tau.qPCR.s <- function(df) {
  acf.qPCR.s = data.frame(acf(df$qPCR_s, pl = F)$acf)
  names(acf.qPCR.s) <- c("acf")
  acf.qPCR.s <- acf.qPCR.s %>% dplyr::mutate(tau = c(1:6))
  
  for(i in 1:dim(acf.qPCR.s)[1]) {
      if(acf.qPCR.s[i, 1] < 0.5) {
        return(acf.qPCR.s[i, 2])
        break
      }
    }
}

tau.qPCR.as <- function(df) {
  acf.qPCR.as = data.frame(acf(df$qPCR_as, pl = F)$acf)
  names(acf.qPCR.as) <- c("acf")
  acf.qPCR.as <- acf.qPCR.as %>% dplyr::mutate(tau = c(1:6))
  
  for(i in 1:dim(acf.qPCR.as)[1]) {
      if (acf.qPCR.as[i, 1] < 0.5) {
        return(acf.qPCR.as[i, 2])
        break
      }
    }
}

tau.dPCR.s <- function(df) {
  acf.dPCR.s = data.frame(acf(df$dPCR_s, pl = F)$acf)
  names(acf.dPCR.s) <- c("acf")
  acf.dPCR.s <- acf.dPCR.s %>% dplyr::mutate(tau = c(1:6))
  
  for(i in 1:dim(acf.dPCR.s)[1]) {
      if (acf.dPCR.s[i, 1] < 0.5) {
        return(acf.dPCR.s[i, 2])
        break
      }
    }
}

tau.dPCR.as <- function(df) {
  acf.dPCR.as = data.frame(acf(df$dPCR_as, pl = F)$acf)
  names(acf.dPCR.as) <- c("acf")
  acf.dPCR.as <- acf.dPCR.as %>% dplyr::mutate(tau = c(1:6))
  
  for(i in 1:dim(acf.dPCR.as)[1]) {
      if (acf.dPCR.as[i, 1] < 0.5) {
        return(acf.dPCR.as[i, 2])
        break
      }
    }
}

noise.space.df <- function(clone, cv) {
  tau.qPCR.s.li.o <- tau.qPCR.s(clone)
  tau.qPCR.as.li.o <- tau.qPCR.s(clone)
  tau.dPCR.s.li.o <- tau.qPCR.s(clone)
  tau.dPCR.as.li.o <- tau.qPCR.s(clone)
  
  df.tau <- data.frame(tau.qPCR.s.li.o, tau.qPCR.as.li.o, tau.dPCR.s.li.o, tau.dPCR.as.li.o)
  names(df.tau) <- c("tau.qPCR.s", "tau.qPCR.as", "tau.dPCR.s", "tau.dPCR.as")
  
  df.full <- dplyr::bind_cols(cv, df.tau)
  return(df.full)
}
