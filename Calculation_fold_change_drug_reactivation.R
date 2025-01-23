##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: Measures of FACS-related parameters.
# Object: Calculation of the fold change for each parameters after sinpro clones were treated by drugs.
##############################################################
# R functions

drug_enrichment_foldchange <- function(df) {
  Mock_1 <- df %>% dplyr::filter(Drug == "Mock_1") %>% dplyr::select(GFP_pos, Mean, CV, Mean.live, CV.live)
  m.value.Mock_1 <- c(mean(Mock_1$GFP_pos), mean(Mock_1$Mean), mean(Mock_1$CV*Mock_1$CV), mean(Mock_1$Mean.live), mean(Mock_1$CV.live*Mock_1$CV.live))
  Mock_2 <- df %>% dplyr::filter(Drug == "Mock_2") %>% dplyr::select(GFP_pos, Mean, CV, Mean.live, CV.live)
  m.value.Mock_2 <- c(mean(Mock_2$GFP_pos), mean(Mock_2$Mean), mean(Mock_2$CV*Mock_2$CV), mean(Mock_2$Mean.live), mean(Mock_2$CV.live*Mock_2$CV.live))
  DMSO <- df %>% dplyr::filter(Drug == "DMSO") %>% dplyr::select(GFP_pos, Mean, CV, Mean.live, CV.live)
  m.value.DMSO <- c(mean(DMSO$GFP_pos), mean(DMSO$Mean), mean(DMSO$CV*DMSO$CV), mean(DMSO$Mean.live), mean(DMSO$CV.live*DMSO$CV.live))
  PMA_I <- df %>% dplyr::filter(Drug == "PMA_I") %>% dplyr::select(GFP_pos, Mean, CV, Mean.live, CV.live)
  m.value.PMA_I <- c(mean(PMA_I$GFP_pos), mean(PMA_I$Mean), mean(PMA_I$CV*PMA_I$CV), mean(PMA_I$Mean.live), mean(PMA_I$CV.live*PMA_I$CV.live))
  PEP <- df %>% dplyr::filter(Drug == "PEP") %>% dplyr::select(GFP_pos, Mean, CV, Mean.live, CV.live)
  m.value.PEP <- c(mean(PEP$GFP_pos), mean(PEP$Mean), mean(PEP$CV*PEP$CV), mean(PEP$Mean.live), mean(PEP$CV.live*PEP$CV.live))
  
  Mock = m.value.Mock_2/m.value.Mock_1
  DMSO = m.value.DMSO/m.value.Mock_2
  PMA_I = m.value.DMSO/m.value.PMA_I
  PEP = m.value.DMSO/m.value.PEP
  
  parameter <- c("delta.GFP_pos", "delta.Mean", "delta.CV", "delta.Mean.live", "delta.CV.live")
  
  df <- t(data.frame(parameter, Mock, DMSO, PMA_I, PEP))
  colnames(df) <- df[1,]
  df <- df[-1,]
  df <- as.data.frame(df)
  df <- rownames_to_column(df, var = "drug")
  
  return(df)
}
