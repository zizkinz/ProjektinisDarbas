bed <- read.table("/Users/maksimcizov/PycharmProjects/ProjektinisDarbas/Genes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
bed <- unique(bed)
bed$V1 <- sub("^chr", "", bed$V1)
bed <- bed[bed$V3 - bed$V2 >101,]


uq <- unique(bed$V4)
dt <- data.frame(V1 = character(), V2 = integer(), V3 = integer(), V4 = character())
for(k in 1:length(uq)){
 l <- bed[bed$V4 == uq[k],]
 # print(list(l$V1[1], min(l$V2), max(l$V3), l$V4[1]))
 print(l$V4[1])
 dt[k,] <-  list(l$V1[1], min(l$V2), max(l$V3), l$V4[1])
}

dt <- dt[dt$V1 != "Un_GL456372v1",]

write.table(dt, file = "gens.bed", sep = "\t", quote = FALSE, row.names = FALSE)
