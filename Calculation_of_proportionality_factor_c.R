#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2024.12.28
##############################################################
# Input: Parameters: qPCR and dPCR measures of sense- and antisense RNA and corresponding (%) GFP-positive cells, population-based CV and the autocorrelation time per sinpro clone.
# Object: Calculation of the proportionality factor (C)
##############################################################

## Mean fluoresence intensity (MFI)
Table_S2 <- read.table("/Users/hchen/Documents/Praca/Projekten/Moi/HIV/Sensitized_HIV_latency_model/df_wichtig/Noise_space/sinpro_clones/Table_S2_noise_space_2025_v3.csv", header = T, stringsAsFactors = F)

group_A <- Table_S2 %>% dplyr::filter(group == "A")
MFI.groupA <- mean(group_A$Mean.live)
group_Am <- Table_S2 %>% dplyr::filter(group == "A-")
MFI.groupAm <- mean(group_Am$Mean.live)
group_S <- Table_S2 %>% dplyr::filter(group == "S")
MFI.groupS <- mean(group_S$Mean.live)

# R functions

cal_parameter_c_pb.live.cv.qPCR <- function(df) {
  c_sRNA <- data.frame((df$MFI)/(df$sRNA.qPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.sRNA.qPCR)))
  names(c_sRNA) <- c("parameter_c")
  c_asRNA <- data.frame((df$MFI)/(df$asRNA.qPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.asRNA.qPCR)))
  names(c_asRNA) <- c("parameter_c")
  c_RNAs_sub <- data.frame(((df$MFI)/(df$sRNA.qPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.sRNA.qPCR)))-((df$MFI)/(df$asRNA.qPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.asRNA.qPCR))))
  names(c_RNAs_sub) <- c("parameter_c")
  c_RNAs <- data.frame((df$MFI)/(df$sRNA.qPCR - df$asRNA.qPCR)*(1/df$live_noise_T1.2/((1/df$live_noise_T1.2 + 1/df$tau.sRNA.qPCR)-(1/df$live_noise_T1.2 + 1/df$tau.asRNA.qPCR))))
  names(c_RNAs) <- c("parameter_c")
  
  df <- bind_rows(c_sRNA, c_asRNA, c_RNAs_sub, c_RNAs)
  df$type <- c(rep("sRNA", 8), rep("asRNA", 8), rep("subtraction", 8), rep("reshaped", 8))
  df$group <- rep(c("S", "S", "A_min", "A_min", "S", "A_min", "A", "S"), 4)
  df$clone <- rep(c("778", "788", "789", "795", "800", "802", "812", "815"), 4)
  df$method <- c(rep("qPCR", 32))
  
  return(df)
}

cal_parameter_c_pb.live.cv.dPCR <- function(df) {
  c_sRNA <- data.frame((df$MFI)/(df$sRNA.dPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.sRNA.dPCR)))
  names(c_sRNA) <- c("parameter_c")
  c_asRNA <- data.frame((df$MFI)/(df$asRNA.dPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.asRNA.dPCR)))
  names(c_asRNA) <- c("parameter_c")
  c_RNAs_sub <- data.frame(((df$MFI)/(df$sRNA.dPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.sRNA.dPCR)))-((df$MFI)/(df$asRNA.qPCR)*(1/df$live_noise_T1.2/(1/df$live_noise_T1.2 + 1/df$tau.asRNA.dPCR))))
  names(c_RNAs_sub) <- c("parameter_c")
  c_RNAs <- data.frame((df$MFI)/(df$sRNA.dPCR - df$asRNA.dPCR)*(1/df$live_noise_T1.2/((1/df$live_noise_T1.2 + 1/df$tau.sRNA.dPCR)-(1/df$live_noise_T1.2 + 1/df$tau.asRNA.dPCR))))
  names(c_RNAs) <- c("parameter_c")
  
  df <- bind_rows(c_sRNA, c_asRNA, c_RNAs_sub, c_RNAs)
  df$type <- c(rep("sRNA", 8), rep("asRNA", 8), rep("subtraction", 8), rep("reshaped", 8))
  df$group <- rep(c("S", "S", "A_min", "A_min", "S", "A_min", "A", "S"), 4)
  df$clone <- rep(c("778", "788", "789", "795", "800", "802", "812", "815"), 4)
  df$method <- c(rep("dPCR", 32))
  
  return(df)
}
