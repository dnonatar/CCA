# Both data must have the same number of samples

microb = read.table("/home/ratanond/Desktop/Thesis/Synthetic/Arghavan/MetaSample2/samples.txt")
microb = as.matrix(microb)
microb_rows = dim(microb)[1]/2    # choose only the first class
microb = microb[1:microb_rows,-1]
#head(t(microb))

rnaseq = read.table("/home/ratanond/Desktop/Thesis/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/easyparams_trn.txt",sep = ",")
rnaseq = as.matrix(rnaseq)
rnaseq_rows = dim(rnaseq)[1]/2  
rnaseq = rnaseq[1:rnaseq_rows,]
rnaseq = rnaseq[,1:150]
#head(t(rnaseq))

library(PMA)
set.seed(1105)
ccaPermImmun = CCA.permute(x = microb, z = rnaseq,
                           typex = "standard", typez = "standard", 
                           nperms = 30, niter = 5, standardize = F)
penXtemp = ccaPermImmun$bestpenaltyx
penZtemp = ccaPermImmun$bestpenaltyz
ccaRsltImmun = CCA(x = microb, z = rnaseq,
                   typex = "standard", typez = "standard",
                   penaltyx = penXtemp, penaltyz = penZtemp,
                   K = 2, niter = 5, standardize = F)
sum(ccaRsltImmun$u != 0)
sum(ccaRsltImmun$v != 0)

ccaScoreU = microb %*% ccaRsltImmun$u
ccaScoreV = rnaseq %*% ccaRsltImmun$v
ccaScores = cbind(ccaScoreU, ccaScoreV)
colnames(ccaScores) = c("U1", "U2", "V1", "V2")
ccaScores = as.data.frame(ccaScores)
ccaScores$type = c(rep("type1", 17), rep("type2", 17))


library(ggplot2)
myCCAPlot = function(x = U1, y = U2, col = V1, shape = type, data = ccaScores,
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
    scale_fill_discrete(name = "Feeding Type",
                        labels = c("type1", "type2")) +
    scale_shape_discrete(name = "Feeding Type",
                         labels = c("type1", "type2"))
  myPlot1
}
myCCAPlot()
myCCAPlot(V1, V2, U1, xyName = "rnaseq", coloName = "microbial")

