### Start to do sparse PCA for Immunology and SEEDlevel2 ###
gaussian = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/params_sb_high_1_ma_trn.txt",sep = ",")
gaussian <- as.matrix(gaussian)
gaussian <- gaussian[1:34,]

set.seed(1105)
spcaImmunCv = SPC.cv(x = gaussian, sumabsvs = seq(1, 5, by = 0.2),
                     nfolds = 10, niter = 5)
spcaImmunRslt = SPC(x = gaussian, sumabsv = spcaImmunCv$bestsumabsv,
                    K = 2, niter = 5)
spcaImmunScores = as.data.frame(spcaImmunRslt$u)
sum(spcaImmunRslt$v != 0)
#rownames(spcaImmunScores) = colnames(immunGeneExpr)
names(spcaImmunScores) = c("Host1", "Host2")
spcaImmunScores$type = c(rep("BF", 34), rep("FF", 0))

## Sparse PCA plot for Immunology gene expression levels:
myPCAPlot = function(x = Host1, y = Host2, shapeCol = type,
                     data = spcaImmunScores, shapeColName = "Feeding Type",
                     xyName = "Host", titleName = "Immunology of Host",
                     textVjust = -1.0, elliLev = 0.6, ...){
  jitterPara = list(...)
  x = deparse(substitute(x))
  y = deparse(substitute(y))
  shapeCol = deparse(substitute(shapeCol))
  if(!"height" %in% names(jitterPara)){
    jitterPara = c(jitterPara, height = 0.01) 
  } else if(!"width" %in% names(jitterPara)){
    jitterPara = c(jitterPara, width = 0.01)
  } 
  spcaPlot = ggplot(data, aes(x = data[,x], y = data[,y],
                              shape = data[,shapeCol], colour = data[,shapeCol])) +
    geom_point(size = 4) +
    scale_shape_discrete(name = shapeColName,
                         labels = c("Breastfeeding", "Formula Feeding")) +
    scale_color_discrete(name = shapeColName,
                         labels = c("Breastfeeding", "Formula Feeding")) +
    geom_text(aes(label = rownames(data)), col = "black",
              vjust = textVjust, size = 4,
              position = do.call("position_jitter", args = jitterPara)) +
    ## The position_jitter will make the values within a group 
    ## a litter bit separate.
    ## On ther other hand, position_dodge will separate the values between groups.
    scale_x_continuous(paste0("First PCA scores from ", xyName)) +
    scale_y_continuous(paste0("Second PCA scores from ", xyName)) +
    labs(title = paste0("Sparse PCA Scores from ", titleName)) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 14, vjust = 1.1, face = "bold")) +
    stat_ellipse(aes(fill = data[, shapeCol]),
                 level = elliLev, geom = "polygon", alpha = 0.2, linetype = 0) +
    scale_fill_discrete(name = shapeColName,
                        labels=c("Breastfeeding","Formula Feeding"))
  spcaPlot
}
myPCAPlot()

## Sparse PCA for Microbial Seedlevl2 

set.seed(1105)
spcaSeedLev2Cv = SPC.cv(x = t(seedLev2TestNorm), sumabsvs = seq(1, 5, by = 0.2),
                        nfolds = 10, niter = 5)
spcaSeedLev2Rslt = SPC(x = t(seedLev2TestNorm), sumabsv = spcaSeedLev2Cv$bestsumabsv,
                       K = 2)
spcaSeedLev2Scores = as.data.frame(spcaSeedLev2Rslt$u)
sum(spcaSeedLev2Rslt$v !=0)
rownames(spcaSeedLev2Scores) = colnames(seedLev2TestNorm)
names(spcaSeedLev2Scores) = c("Host1", "Host2")
spcaSeedLev2Scores$type = c(rep("BF", 6), rep("FF", 6))

## Sparse PCA plot for Microbial Seedlevl2
myPCAPlot(data = spcaSeedLev2Scores, xyName = "SeedLevel2",
          titleName = "Microbial Seedlevel2")