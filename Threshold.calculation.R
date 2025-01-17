##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: Measures of qPCR and dPCR from sense and antisense RNAs.
# Object: Determination of the threshold of sense and antisense RNA transcriptions across sinpro clones in different groups.
##############################################################
# R functions

Box.features.calculate <- function(df) {
  df.qPCR.s <- data.frame(print(boxplot.stats(df$qPCR_s)$stats))
  names(df.qPCR.s) <- c("measure")
  df.qPCR.s <- dplyr::mutate(df.qPCR.s, method = "qPCR", RNA = "sRNAs", box.fi = c("Min", "Q1", "Q2", "Q3", "Max"))
  
  df.qPCR.as <- data.frame(print(boxplot.stats(df$qPCR_as)$stats))
  names(df.qPCR.as) <- c("measure")
  df.qPCR.as <- dplyr::mutate(df.qPCR.as, method = "qPCR", RNA = "asRNAs", box.fi = c("Min", "Q1", "Q2", "Q3", "Max"))
  
  df.dPCR.s <- data.frame(print(boxplot.stats(df$dPCR_s)$stats))
  names(df.dPCR.s) <- c("measure")
  df.dPCR.s <- dplyr::mutate(df.dPCR.s, method = "dPCR", RNA = "sRNAs", box.fi = c("Min", "Q1", "Q2", "Q3", "Max"))
  
  df.dPCR.as <- data.frame(print(boxplot.stats(df$dPCR_as)$stats))
  names(df.dPCR.as) <- c("measure")
  df.dPCR.as <- dplyr::mutate(df.dPCR.as, method = "dPCR", RNA = "asRNAs", box.fi = c("Min", "Q1", "Q2", "Q3", "Max"))
  
  df.pool <- dplyr::bind_rows(df.qPCR.s, df.qPCR.as, df.dPCR.s, df.dPCR.as)
  
  return(df.pool)
}
