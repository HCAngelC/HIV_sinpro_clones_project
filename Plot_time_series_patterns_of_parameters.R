##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: Master table with FACS-related variables and measures from sRNAs and asRNAs.
# Object: Plot the time-series pattern of FACS-related variables and transcription of sRNAs and asRNAs over time.
##############################################################
# R functions

Plot_individual_phenotype_overTime <- function(df) {
  p_gfp <- ggplot(df, aes(x = Date, y = GFP_pos, group = 1))+geom_line(color = "#66CC00")+geom_point(color = "#66CC00" )+theme_bw()+xlab("Time")+ylab("(%) GFP(+)")+ylim(0, 100)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_mean <- ggplot(df, aes(x = Date, y = log2(Mean.live), group = 1))+geom_line(color = "#996600")+geom_point(color = "#996600" )+theme_bw()+xlab("Time")+ylab("Mean.live (log2)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_CV_square <- ggplot(df, aes(x = Date, y = log2(CV.live*CV.live), group = 1))+geom_line(color = "orange")+geom_point(color = "orange" )+theme_bw()+xlab("Time")+ylab("CV.live square (log2)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_qPCR_s <- ggplot(df, aes(x = Date, y = qPCR_s, group = 1))+geom_line(color = "#CC0066")+geom_point(color = "#CC0066" )+theme_bw()+xlab("Time")+ylab("Sense RNA (qPCR)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_qPCR_as <- ggplot(df, aes(x = Date, y = qPCR_as, group = 1))+geom_line(color = "#FF66FF")+geom_point(color = "#FF66FF" )+theme_bw()+xlab("Time")+ylab("Antisense RNA (qPCR)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_dPCR_s <- ggplot(df, aes(x = Date, y = dPCR_s, group = 1))+geom_line(color = "#0000FF")+geom_point(color = "#0000FF" )+theme_bw()+xlab("Time")+ylab("Sense RNA (dPCR)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  p_dPCR_as <- ggplot(df, aes(x = Date, y = dPCR_as, group = 1))+geom_line(color = "#3399FF")+geom_point(color = "#3399FF" )+theme_bw()+xlab("Time")+ylab("Antisense RNA (dPCR)")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=8, colour = "black",angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 8, colour = "black"))
  
  f <- ggarrange(p_gfp, p_mean, p_CV_square, p_qPCR_s, p_qPCR_as, p_dPCR_s, p_dPCR_as, ncol = 1, nrow = 7)
  print(f)
}
