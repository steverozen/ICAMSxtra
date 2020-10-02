# Source this file from ICAMSxtra top level directory.

cat(getwd(), "\n")

files <- list.files(path = "data-raw/vcfs/HepG2_Cis/", full.names = TRUE)
vcf.names <- gsub(".vcf", "", tools::file_path_sans_ext(basename(files)))
HepG2.catalogs <- 
  ICAMS::StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                     region = "genome",
                                     names.of.VCFs = vcf.names)
HepG2.SBS96 <- HepG2.catalogs$catSBS96
WriteCatalog(HepG2.SBS96, file = "tests/testthat/testdata/HepG2_Cis.spectra.csv")

files1 <- list.files(path = "data-raw/vcfs/MCF10A_Cis/", full.names = TRUE)
vcf.names1 <- gsub(".vcf", "", tools::file_path_sans_ext(basename(files1)))
MCF10A.catalogs <- 
  ICAMS::StrelkaSBSVCFFilesToCatalog(files1, ref.genome = "hg19", 
                                     region = "genome",
                                     names.of.VCFs = vcf.names1)
MCF10A.SBS96 <- MCF10A.catalogs$catSBS96
WriteCatalog(MCF10A.SBS96, file = "tests/testthat/testdata/MCF10A_Cis.spectra.csv")
