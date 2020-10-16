GRCh37.proportions <- list()
# total size of genome for GRCh37
# calculated from loading the BSgenome object for GRCh37 and summing the length of
# all chromosomes (including X and Y)
GRCh37.proportions$total.bp <- 3095677412

# size of genic portion of genome
GRCh37.proportions$coding.bp <- 
  Reduce("+", ICAMS::trans.ranges.GRCh37$end - ICAMS::trans.ranges.GRCh37$start + 1)

# size of the intergenic portion of genome
GRCh37.proportions$noncoding.bp <- 
  GRCh37.proportions$total.bp - GRCh37.proportions$coding.bp

# proportion of genic portion of genome
GRCh37.proportions$prop.coding <- 
  GRCh37.proportions$coding.bp / GRCh37.proportions$total.bp

# proportion of intergenic portion of genome
GRCh37.proportions$prop.noncoding <- 
  GRCh37.proportions$noncoding.bp / GRCh37.proportions$total.bp


GRCh38.proportions <- list()
# total size of genome for GRCh38
# example of how to calculate total.bp for hg38
# hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# tbl <- BSgenome::seqinfo(hg38)
# total.bp <- sum(tbl@seqlengths[1:24])
GRCh38.proportions$total.bp <- 3088269832

# size of genic portion of genome
GRCh38.proportions$coding.bp <- Reduce("+", ICAMS::trans.ranges.GRCh38$end - ICAMS::trans.ranges.GRCh38$start + 1)

# size of the intergenic portion of genome
GRCh38.proportions$noncoding.bp <- 
  GRCh38.proportions$total.bp - GRCh38.proportions$coding.bp

# proportion of genic portion of genome
GRCh38.proportions$prop.coding <- 
  GRCh38.proportions$coding.bp / GRCh38.proportions$total.bp

# proportion of intergenic portion of genome
GRCh38.proportions$prop.noncoding <- 
  GRCh38.proportions$noncoding.bp / GRCh38.proportions$total.bp
