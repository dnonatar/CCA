# Both data must have the same number of samples

gaussian <- read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/myparams_ma_trn.txt",sep = ",")
gaussian <- as.matrix(gaussian)
gaussian_rows <- dim(gaussian)[1]/2  
gaussian <- gaussian[1:gaussian_rows,1:3] # try with 3 genes

# choose OTU that has the highest mean
# y is y_max without the zeros
OTU_max <- which.max(apply(microb,2,mean))
OTU_max_2 <- which.max(apply(microb[,-142],2,mean))
microb <- microb[,c(142,147)]        # try with 2 OTUs
y_max <- microb[,which.max(apply(microb,2,mean))]
y <- y_max[y_max>0]

# choose genes from the same block (block size = 5)
gaussian_block <- gaussian[,1:2]
x_all_1 <- gaussian_block[,1]
x_all_2 <- gaussian_block[,2]
x1 <- x_all_1[which(y>0)]
x2 <- x_all_2[which(y>0)]

model <- lm(y~x1+x2)
summary(model)
intercept <- model$coefficients[1]
slope1 <- model$coefficients[2]
slope2 <- model$coefficients[3]

#xx <- -1:1
#yy <- slope*xx+intercept
#plot(x,y, xlab = "gaussian", ylab = "microb")
#lines(xx,yy,type = "l")

#x_adjust <- x[order(x)]
#y_adjust <- slope1*x1+slope2*x2+intercept
#plot(x_adjust,y_adjust, xlab = "gaussian", ylab = "microb")
#lines(xx,yy,type = "l")

y_new <- y_max

for (i in 1:length(y_new)) {
  if (y_new[i]>0) {
    y_new[i]<- slope1*x_all_1[i]+slope2*x_all_2[i]+intercept
  }
}

microb_new <- microb
microb_new[,1] <- y_new


library(PMA)
set.seed(1105)
ccaPerm = CCA.permute(x = microb_new, z = rnaseq,
                           typex = "standard", typez = "standard", 
                           nperms = 30, niter = 5, standardize = T)
penXtemp = ccaPerm$bestpenaltyx
penZtemp = ccaPerm$bestpenaltyz
ccaRslt = CCA(x = microb_new, z = rnaseq,
                   typex = "standard", typez = "standard",
                   penaltyx = penXtemp, penaltyz = penZtemp,
                   K = 2, niter = 5, standardize = T)
sum(ccaRslt$u != 0)
sum(ccaRslt$v != 0)

ccaScoreU = microb_new %*% ccaRslt$u
ccaScoreV = rnaseq %*% ccaRslt$v
ccaScores = cbind(ccaScoreU, ccaScoreV)
colnames(ccaScores) = c("U1", "U2", "V1", "V2")
ccaScores = as.data.frame(ccaScores)
#number of each type should be flexible
ccaScores$type = c(rep("class 1", 17), rep("class 2", 17))

# find positions of the top 5 absolute values for U1
match(sort(abs(ccaRslt$u),decreasing = T)[1:5],abs(ccaRslt$u))
match(sort(abs(ccaRslt_old$u),decreasing = T)[1:5],abs(ccaRslt_old$u))

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
    scale_fill_discrete(name = "Class",
                        labels = c("class 1", "class 2")) +
    scale_shape_discrete(name = "Class",
                         labels = c("class 1", "class 2"))
  myPlot1
}
myCCAPlot()
myCCAPlot(V1, V2, U1, xyName = "rnaseq", coloName = "microbial")


