# install.packages("devtools")
# devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
# 
# BiocManager::install(c("rtracklayer","GenomicRanges"))
# 
# LGG_MAF <- read.delim(
#   file = "../vignettes/TCGA.LGG.mutect.2c0dab06-7af9-4380-bc83-13148298da19.DR-7.0.somatic.maf",
#   header = T,skip = 5,stringsAsFactors = F)
# 
# BiocManager::install("GenomicRanges")
# library(GenomicRanges)
# library(rtracklayer)
# 
# biocLite("BSgenome")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# install.packages("deconstructSigs",repos = "https://cloud.r-project.org")