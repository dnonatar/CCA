# Both data must have the same number of samples

gaussian <- read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/myparams_ma_trn.txt",sep = ",")
gaussian <- as.matrix(gaussian)
gaussian_rows <- dim(gaussian)[1]/2  
gaussian <- gaussian[1:gaussian_rows,1:3] # try with 3 genes

# choose OTUs that have the smallest number of zeros
f <- function(m){return(sum(m!= 0))}
order(apply(microb,2,f),decreasing = T)[1:10]  
microb <- microb[,c(417,541,42)]        # try with 2 OTUs
#y_max <- microb[,which.max(apply(microb,2,mean))]
#y <- y_max[y_max>0]   # y is y_max without the zeros
y <- microb[which(microb[,1]>0),1]

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

#y_new <- y_max
y_new <- microb[,1]

for (i in 1:length(y_new)) {
  #if (y_new[i]>0) {
    y_new[i]<- slope1*x_all_1[i]+slope2*x_all_2[i]+intercept
 # }
}

microb_new <- microb
microb_new[,1] <- y_new


