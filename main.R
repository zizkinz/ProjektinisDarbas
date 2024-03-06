library(CoverageView)
library(MASS)

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

fn <- function (x) {
  x <- na.omit(x)
  return (-mean(x))
}


ngenes <- 8


# -17 -11 -10
# filelist <- c("test_data/DP5.bam","test_data/DP12.bam", "test_data/DP13.bam", "test_data/DP14.bam", "test_data/DP15.bam", "test_data/DP16.bam", "test_data/DP18.bam", "test_data/DP19.bam", "test_data/DP20.bam", "test_data/DP21.bam", "test_data/DP22.bam", "test_data/DP23.bam", "test_data/DP24.bam", "test_data/DP25.bam")

nsamples <- 25
ks_mx <- matrix(ncol = ngenes, nrow = nsamples)
nms <- as.data.frame(read.table("test_data/Genes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))$V4
colnames(ks_mx) <- nms
for (j in 1:nsamples){
  trm<-CoverageBamFile(paste("test_data/DP",j,".bam",sep = ''))
  # trm<-CoverageBamFile("test_data/DP23_sequence.bam")
  # trm<-CoverageBamFile("test_data/DP24_sequence.bam")

  n_wind <- 100
  covmx<-cov.matrix(trm, coordfile= "test_data/Genes.bed", bin_width=1, no_windows = n_wind)

  # write.profile(DP23_Glul_covmx,outfile="DP23_sequence_Glu.txt")
  # draw.profile(covmx,ylab="avg coverage",outfile=paste(j, ".png", sep = ''))
  # draw.heatmap(covmx,outfile="DP23_sequence_Glu_heatmap.png")

  # plot(1:nrow(covmx),cumsum(covmx[,6]/sum(covmx[,6])), type = "l")
  # lines(1:nrow(covmx), punif(1:nrow(covmx),0,100), type = "l")

  ks_mx[j,] <- ks(hyp_dist, covmx)
}

ks_norm <- ks_mx - median(ks_mx,na.rm = TRUE)
mrin <- apply(ks_norm, MARGIN = 1, FUN = fn)
print(mrin)

# fitting mrins
hist(mrin)

fitdistr(mrin,densfun = "normal")$sd

#       mean         sd
# 0.01935427 0.01368553

qqnorm(mrin, pch = 16)
qqline(mrin)

# error <- qnorm(0.995) * as.numeric(l[2])/sqrt(n_wind)

conf_int <- c(as.numeric(l[1]) - error, as.numeric(l[1]) + error)

plot(ecdf(mrin))
abline(v=l[1], col = "red")
abline(v=l[1] - 2 * l[2], col = "blue")
abline(v=l[1] + 2 * l[2], col = "blue")


#   dat <- read.table("test_data/DP23_sequence_Glul.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
# hist((dat$V2), breaks = 100, freq = F)
# hist(runif(nrow(dat), min = min(dat$V2), max = max(dat$V2)))
# ks.test(dat$V2, 'punif',min(dat$V2), max(dat$V2))


