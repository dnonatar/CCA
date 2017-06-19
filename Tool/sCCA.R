# Both data must have the same number of samples

microb = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Arghavan/MetaSample2/samples.txt")
microb = as.matrix(microb)
microb <- microb[,-1]
microb_rows = dim(microb)[1]/2   # divide by two to choose only the first class
microb = microb[1:microb_rows,] 
sd <- apply(microb,2,sd)
microb <- microb[,which(sd!=0)]  # choose only OTUs with non-zero standard deviation
#head(t(microb))

rnaseq = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/params_sb_high_1_trn.txt",sep = ",")
rnaseq = as.matrix(rnaseq)
rnaseq_rows = dim(rnaseq)[1]/2  
rnaseq = rnaseq[1:rnaseq_rows,1:600]  
#rnaseq = rnaseq[,1:150]
#head(t(rnaseq))

gaussian = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/params_sb_high_1_ma_trn.txt",sep = ",")
gaussian = as.matrix(gaussian)
gaussian_rows = dim(rnaseq)[1]/2  
gaussian = rnaseq[1:rnaseq_rows,]
gaussian = rnaseq[,1:150]

library(PMA)
set.seed(1105)
ccaPerm_old = CCA.permute(x = microb, z = rnaseq,
                           typex = "standard", typez = "standard", 
                           nperms = 30, niter = 5, standardize = T)
penXtemp = ccaPerm_old$bestpenaltyx
penZtemp = ccaPerm_old$bestpenaltyz
ccaRslt_old = CCA(x = microb, z = rnaseq,
                   typex = "standard", typez = "standard",
                   penaltyx = penXtemp, penaltyz = penZtemp,
                   K = 2, niter = 5, standardize = T)
sum(ccaRslt_old$u != 0)
sum(ccaRslt_old$v != 0)

ccaScoreU_old = microb %*% ccaRslt_old$u
ccaScoreV_old = rnaseq %*% ccaRslt_old$v
ccaScores_old = cbind(ccaScoreU_old, ccaScoreV_old)
colnames(ccaScores_old) = c("U1", "U2", "V1", "V2")
ccaScores_old = as.data.frame(ccaScores_old)
#number of each type should be flexible
ccaScores_old$type = c(rep("class 1", 17), rep("class 2", 17))


library(ggplot2)
myCCAPlot = function(x = U1, y = U2, col = V1, shape = type, data = ccaScores_old,
                     xyName = "microbial", coloName = "rnaseq",
                     textVjust = -1.0, elliLev = 0.6, ...){
  jitterPara = list(...)
  if(!"height" %in% names(jitterPara)){
    jitterPara = c(jitterPara, height = 0.01) 
  } else if(!"width" %in% names(jitterPara)){
    jitterPara = c(jitterPara, width = 0.01)
  }
  x = deparse(substitute(x))
  y = deparse(substitute(y))
  col = deparse(substitute(col))
  shape = deparse(substitute(shape))
  myPlot1 = ggplot(data, aes(x = data[,x], y = data[,y],
                             col = data[,col], shape = data[,shape])) +
    geom_point(size = 4) +
    scale_color_continuous(name = paste0("First Component \nScores of ",
                                         coloName),
                           low = "blue", high = "red") +
    geom_text(aes(label = rownames(data)),
              col = "black", size = 5, vjust = textVjust,
              position = do.call("position_jitter", args = jitterPara)) +
    ## The position_jitter will make the values within a group 
    ## a litter bit separate.
    ## On ther other hand, position_dodge will separate the values between groups.
    scale_x_continuous(paste0("First Component Scores of ",
                              xyName)) +
    scale_y_continuous(paste0("Second Component Scores of ",
                              xyName)) +
    labs(title = paste0("Sparse CCA Scores for ", xyName, " as Base")) +
    theme(legend.title = element_text(size = 12),
          plot.title = element_text(size = 16, vjust = 2.0, face = "bold"),
          legend.text = element_text(size = 10)) +
    stat_ellipse(aes(fill = data[,shape]), level = elliLev, alpha = 0.2,
                 geom = "polygon", linetype = 2) +
    scale_fill_discrete(name = "Class",
                        labels = c("class 1", "class 2")) +
    scale_shape_discrete(name = "Class",
                         labels = c("class 1", "class 2"))
  myPlot1
}


## Write output files (1 csv + 2 pdf)

pdf(file.path("/home/ratanond/Desktop/Masters_Project/CCA/Tool/ccaresults",paste(randname,"Scores1.pdf", sep="")),bg="transparent")
myCCAPlot()
dev.off()
pdf(file.path("/home/ratanond/Desktop/Masters_Project/CCA/Tool/ccaresults",paste(randname,"Scores2.pdf", sep="")),bg="transparent")
myCCAPlot(V1, V2, U1, xyName = "Immunology", coloName = "Seedlevel2")
dev.off()

outfile = file.path("/home/ratanond/Desktop/Masters_Project/CCA/Tool/ccaresults",  
                    paste(randname,"-", "ccaImmunScores", ".csv",sep=""))
#write.csv(x = ccaScores_old, file = outfile)
