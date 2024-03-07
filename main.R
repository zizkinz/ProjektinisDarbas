library(CoverageView)
library(MASS)
library(plot.matrix)
library(RColorBrewer)

hyp_dist <- function(x){
  return(punif(x,0,n_wind))
}

ks <- function(hyp_dist, covmx){
  ks_d <- numeric(ncol(covmx))
  for (p in 1:ncol(covmx)){
  emp_dist <- function(x){
    cv <- cumsum(covmx[,p]/sum(covmx[,p]))
    if (x < 0){
      return(0)
    }
    else if (x > n_wind){
    return(1)
}
    else {
      return(cv[as.integer(x)])
    }

}
  d <- numeric(n_wind)
  for (k in (1:n_wind)){
    d[k] <- abs(emp_dist(k) - hyp_dist(k))
  }
  ks_d[p] <- max(d)
}
  return(ks_d)
}

ngenes <- 8
nsamples <- 25
n_wind <- 100

ks_mx <- matrix(ncol = ngenes, nrow = nsamples)
nms <- as.data.frame(read.table("test_data/Genes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))$V4
# k <- as.data.frame(read.table("test_data/AllGenes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(ks_mx) <- nms
for (j in 1:nsamples){
  trm<-CoverageBamFile(paste("test_data/DP",j,".bam",sep = ''))
  # trm<-CoverageBamFile("test_data/DP23_sequence.bam")
  # trm<-CoverageBamFile("test_data/DP24_sequence.bam")


  covmx<-cov.matrix(trm, coordfile= "test_data/Genes.bed", bin_width=1, no_windows = n_wind)

  # write.profile(DP23_Glul_covmx,outfile="DP23_sequence_Glu.txt")
  # draw.profile(covmx,ylab="avg coverage",outfile=paste(j, ".png", sep = ''))
  # draw.heatmap(covmx,outfile="DP23_sequence_Glu_heatmap.png")
  # barplot((covmx[,5]))

  # plot(1:nrow(covmx),cumsum(covmx[,6]/sum(covmx[,6])), type = "l")
  # lines(1:nrow(covmx), punif(1:nrow(covmx),0,100), type = "l")

  ks_mx[j,] <- ks(hyp_dist, covmx)
}

ks_norm <- ks_mx - apply(ks_mx, MARGIN = 2, FUN = function (x) rep(median(x, na.rm = TRUE),times = length(x)))

mrin <- apply(ks_norm, MARGIN = 1, FUN = function (x) -mean(na.omit(x)))
options(scipen=999)
print(mrin)

# mKS matrix
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(ks_norm,asp = FALSE, breaks = c(-1,-0.2,-0.05,-0.02,-0.01,0,0.01,0.02,0.05,0.2,1), col=brewer.pal(10,"RdBu"), na.col="white")

# fitting mrins

l <- fitdistr(mrin,densfun = "normal")$sd
l

#       mean         sd
# 0.01643934 0.01162437

# histogram
hist(mrin,breaks = 10)
mean(mrin)
lines(seq(from = -0.25, to= 0.1, by = 0.001),dnorm(seq(from = -0.25, to= 0.1, by = 0.001),mean = l[1], sd = l[2]),col="red")
abline(v=l[1], col = "blue")

# Q-Q grafikas
qqnorm(mrin)
qqline(mrin)


error <- qnorm(0.975,mean = l[1], sd = l[2])
conf_int <- c(as.numeric(l[1]) - error, as.numeric(l[1]) + error)
which(mrin < conf_int[1])
# 1  2  3  4  6  8  9 10 11 12

shapiro.test(mrin)
#       mean         sd
# 0.01643934 0.01162437

plot(ecdf(mrin))
abline(v=l[1], col = "red")
abline(v=conf_int[1], col = "blue")
abline(v=conf_int[2], col = "blue")


