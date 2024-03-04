library(Rsamtools)
filename <-  'test_data/DP23_sequence.bam'
utils:::format.object_size(file.info(filename)$size, "auto")

bamFile <- BamFile(filename)
seqinfo(bamFile)

df <- as.data.frame(seqinfo(bamFile))
head(df)
#
bamFile
seqlevels(bamFile)
dd
genome(bamFile)

BamViews(bamFile)

which <- GRanges("9:107987096-107993120")
param <- ScanBamParam(which=which, what=scanBamWhat())
sq <- scanBam(bamFile, param=param)

seqlevels(bamFile)


ss <- sq$`9:107987096-107993120`$seq
library(stringr)
z <- str_split_1(as.character(ss[1]), pattern = '')
barplot(table(z))

