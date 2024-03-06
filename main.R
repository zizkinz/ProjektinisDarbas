library(CoverageView)
trm<-CoverageBamFile("test_data/DP23_sequence.bam")
# trm<-CoverageBamFile("test_data/DP24_sequence.bam")

n_wind <- 100
covmx<-cov.matrix(trm, coordfile= "test_data/Genes.bed", bin_width=1, no_windows = n_wind)
# write.profile(DP23_Glul_covmx,outfile="DP23_sequence_Glu.txt")
draw.profile(covmx,ylab="avg coverage",outfile="DP23_sequence_Coverage.png")
# draw.heatmap(covmx,outfile="DP23_sequence_Glu_heatmap.png")

plot(1:nrow(covmx),cumsum(covmx[,2]/sum(covmx[,2])), type = "l")
lines(1:nrow(covmx), punif(1:nrow(covmx),0,100), type = "l")

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


# Reikia atiminet mediana per visus meginius
ks_norm <- ks(hyp_dist, covmx) - median(ks(hyp_dist, covmx))

mrin <- -mean(ks_norm)



# # fitting mrins
# hist(ks_norm)
#
# library(MASS)
#
# l <- fitdistr(ks_norm,densfun = "normal")$sd
#  #      mean          sd
#  #  0.05126011   0.20428576
#  # (0.07222592) (0.05107144)
#
#
#
# qqnorm(ks_norm, pch = 16)
# qqline(ks_norm)
# ks_norm
#
# plot(ecdf(ks_norm))
#
# error <- qnorm(0.995)*as.numeric(l[2])/sqrt(n_wind)



#   dat <- read.table("test_data/DP23_sequence_Glul.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
# hist((dat$V2), breaks = 100, freq = F)
# hist(runif(nrow(dat), min = min(dat$V2), max = max(dat$V2)))
# ks.test(dat$V2, 'punif',min(dat$V2), max(dat$V2))


