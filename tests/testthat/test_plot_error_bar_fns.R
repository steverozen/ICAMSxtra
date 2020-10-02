context("Plot error bar functions")

test_that("MeanOfSpectraAsSig function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  MCF10A.catSBS96 <- ICAMS::ReadCatalog(file = "testdata/MCF10A_Cis.spectra.csv",
                                        ref.genome = "hg19", region = "genome")
  retval <- MeanOfSpectraAsSig(spectra = MCF10A.catSBS96)
  retval2 <- MeanOfSpectraAsSig(spectra = MCF10A.catSBS96, mean.weighted = FALSE)
  
  ICAMS::PlotCatalog(retval$mean.sig)
  ICAMS::PlotCatalog(retval2$mean.sig)
  graphics.off()
  expect_equal(retval$constituent.sigs, retval2$constituent.sigs)
})

test_that("PlotSpectraAsSigsWithUncertainty function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  MCF10A.catSBS96 <- ICAMS::ReadCatalog(file = "testdata/MCF10A_Cis.spectra.csv",
                                        ref.genome = "hg19", region = "genome")
  retval <- PlotSpectraAsSigsWithUncertainty(spectra = MCF10A.catSBS96)
  retval2 <- PlotSpectraAsSigsWithUncertainty(spectra = MCF10A.catSBS96, 
                                              mean.weighted = FALSE)
  graphics.off()
  expect_equal(retval$constituent.sigs, retval2$constituent.sigs)
  unlink("Rplots.pdf")
})