# library(Rsamtools)
# filename <-  'test_data/DP23_sequence.bam'
# utils:::format.object_size(file.info(filename)$size, "auto")
#
# bamFile <- BamFile(filename)
# seqinfo(bamFile)
#
# df <- as.data.frame(seqinfo(bamFile))
# head(df)
# #
# bamFile
# seqlevels(bamFile)
#
# genome(bamFile)
#
# # BamViews(bamFile)
#
# which <- GRanges("9:107987096-107993120")
# param <- ScanBamParam(which=which, what=scanBamWhat())
# sq <- scanBam(bamFile, param=param)
#
# ss <- sq$`9:107987096-107993120`$seq
# library(stringr)
# z <- str_split_1(as.character(ss), pattern = '')
# barplot(table(z)


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("CoverageView")
library(CoverageView)
trm<-CoverageBamFile("test_data/DP23_sequence.bam")

DP23_Glul_covmx<-cov.matrix(trm,coordfile="test_data/DP23_sequence_Glul.bed",extend=1000,bin_width=10)
DP23_Glul_covmx[1:3,1:5]

write.profile(DP23_Glul_covmx,outfile="DP23_sequence_Glul.txt")
draw.profile(DP23_Glul_covmx,ylab="avg coverage",outfile="DP23_sequence_Glul.png")
draw.heatmap(DP23_Glul_covmx,outfile="DP23_sequence_Glul_heatmap.png")

