library(CoverageView)
library(MASS)
library(RColorBrewer)

hyp_dist <- function(x){
  return(punif(x,0,n_wind))
}
emp_dist <- function(x, cv){
    if (x < 1){
      return(as.integer(0))
    }
    else if (x > n_wind){
    return(as.integer(1))
    }
    else {
      return(cv[as.integer(x)])
    }
}

ks <- function(hyp_dist, covmx){
  ks_d <- numeric(ncol(covmx))
  for (p in 1:ncol(covmx)){
      if(sum(covmx[,p]) == 0){
        ks_d[p] <- NaN
        next
      }
      cv <- cumsum(covmx[,p]/sum(covmx[,p]))
      d_plus <- numeric(n_wind + 1)
      d_minus <- numeric(n_wind + 1)
      for (k in (1:n_wind)){
        d_plus[k] <- abs(emp_dist(k, cv) - hyp_dist(k))
      }
          for (k in (1:n_wind)){
        d_minus[k] <- abs(emp_dist(k-1, cv) - hyp_dist(k))
      }
      ks_d[p] <- max(c(d_plus,d_minus))
}
  return(ks_d)
}




coverages <- list.files("Cov_matrices/")
new_f <- list.files("Data/")
dat <- new_f[grepl('.*\\.bam$',new_f)]

ngenes <- 46577
ncov <- length(coverages)
nsamples <- length(dat)
n_wind <- 100

ks_mx <- matrix(ncol = ngenes, nrow = ncov + nsamples)
nms <- as.data.frame(read.table("test_data/Genes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))$V4
# k <- as.data.frame(read.table("test_data/AllGenes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(ks_mx) <- nms
rownames(ks_mx) <- c(coverages,dat)
ks_mx <- data.frame(ks_mx)


for (k in 1:ncov){
  # k <- 16
    start_time <- Sys.time()
  print(paste("Loading ", coverages[k]))
  covmx <- read.table(paste("Cov_matrices/",coverages[k],sep=''),header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
  ks_mx[k,] <- ks(hyp_dist, covmx)
  end_time <- Sys.time()
  print(end_time - start_time)
}

# timesss <- numeric(30)
for (k in 1:nsamples){
  # j <- 1
  # start_time <- Sys.time()
  # trm<-CoverageBamFile(paste("test_data/DP",j,".bam",sep = ''))
  # trm<-CoverageBamFile("test_data/DP23_sequence.bam")
  trm<-CoverageBamFile(paste("Data/",dat[k], sep=''))

  covmx<-cov.matrix(trm, coordfile= "test_data/Genes.bed", bin_width=1, no_windows = n_wind)
  # end_time <- Sys.time()
  # timesss[j] <- end_time - start_time
  colnames(covmx) <- nms
  write.table(covmx, file = paste("Cov_matrices/", sub(".bam$", "", dat)[k], "_cov.bed", sep = ''), sep = "\t", quote = FALSE, row.names = FALSE)
  # write.profile(covmx,outfile=paste(j,"_cov.txt"))
  # draw.profile(covmx,ylab="avg coverage",outfile=paste("Graphs/",j, ".png", sep = ''))
  # draw.heatmap(covmx,outfile="DP23_sequence_Glu_heatmap.png")
  # barplot((covmx[,5]))

  # plot(1:nrow(covmx),cumsum(covmx[,6]/sum(covmx[,6])), type = "l")
  # lines(1:nrow(covmx), punif(1:nrow(covmx),0,100), type = "l")

  ks_mx[ncov + k,] <- ks(hyp_dist, covmx)
}

z <- ks(hyp_dist, covmx)

ks_norm <- ks_mx - apply(ks_mx, MARGIN = 2, FUN = function (x) rep(median(x, na.rm = TRUE),times = length(x)))

mrin <- apply(ks_norm, MARGIN = 1, FUN = function (x) -mean(na.omit(x)))
# mrin <- apply(ks_mx, MARGIN = 1, FUN = function (x) -mean(na.omit(x)))
options(scipen=999)
print(mrin)


ks_stat <- 1.36/sqrt(101)
tin <- apply(ks_mx > ks_stat, MARGIN = 1, FUN = function (x) mean(x, na.rm = TRUE))

# mKS matrix
par(mar=c(5.1, 4.1, 4.1, 4.1))
library(plot.matrix)
plot(ks_norm,asp = FALSE, breaks = c(-1,-0.2,-0.05,-0.02,-0.01,0,0.01,0.02,0.05,0.2,1), col=brewer.pal(10,"RdBu"), na.col="white")
# detach("package:plot.matrix", unload=TRUE)
# unloadNamespace("plot")

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
#  1  2  3  4  5  6  7  8  9 10 11 12 13 16 17 18 19

shapiro.test(mrin)
#       mean         sd
# 0.01643934 0.01162437

plot(ecdf(mrin))
abline(v=l[1], col = "red")
abline(v=conf_int[1], col = "blue")
abline(v=conf_int[2], col = "blue")

png(filename="name.png")
plot(tin,mrin)
dev.off()