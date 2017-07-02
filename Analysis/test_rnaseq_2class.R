rnaseq = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/tst/params_sb_high_1_tst.txt",sep = ",")
rnaseq = as.matrix(rnaseq)
rnaseq_0 <- rnaseq[1:50,]
rnaseq_1 <- rnaseq[501:550,]


library(PMA)
set.seed(1105)
ccaPerm_old = CCA.permute(x = rnaseq_0, z = rnaseq_1,
                          typex = "standard", typez = "standard", 
                          nperms = 30, niter = 5, standardize = T)
penXtemp = ccaPerm_old$bestpenaltyx
penZtemp = ccaPerm_old$bestpenaltyz
ccaRslt_old = CCA(x = rnaseq_0, z = rnaseq_1,
                  typex = "standard", typez = "standard",
                  penaltyx = penXtemp, penaltyz = penZtemp,
                  K = 2, niter = 5, standardize = T)
sum(ccaRslt_old$u != 0)
sum(ccaRslt_old$v != 0)

ccaScoreU_old = rnaseq_0 %*% ccaRslt_old$u
ccaScoreV_old = rnaseq_1 %*% ccaRslt_old$v
ccaScores_old = cbind(ccaScoreU_old, ccaScoreV_old)
colnames(ccaScores_old) = c("U1", "U2", "V1", "V2")
ccaScores_old = as.data.frame(ccaScores_old)
#number of each type should be flexible
#ccaScores_old$type = c(rep("class 1", 17), rep("class 2", 17))


library(ggplot2)
myCCAPlot = function(x = U1, y = U2, col = V1, shape = type, data = ccaScores_old,
                     xyName = "rnaseq_0", coloName = "rnaseq_1",
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
myCCAPlot()
myCCAPlot(V1, V2, U1, xyName = "rnaseq_1", coloName = "rnaseq_0")
