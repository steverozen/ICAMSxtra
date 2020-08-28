context("PlotCatalogToPdf.ID115Catalog")

test_that("PlotID115CatalogToPdf function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  test.files <- c("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.vcf",
                  "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s2.vcf")
  list.of.ID.vcfs <- ReadStrelkaIDVCFs(test.files)
  list.of.catalogs <- VCFsToID115Catalogs(list.of.ID.vcfs, ref.genome = "hg19",
                                         region = "genome")
  catalog.counts <- list.of.catalogs$catalog
  colnames(catalog.counts) <- paste0("Strelka.ID.GRCh37.", 1 : 2)
  out <- PlotID115CatalogToPdf(catalog.counts,
                               file = file.path(tempdir(), "PlotCatID115.counts.test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  #test VCFsToID115CatalogsAndPlotToPdf
  out <- VCFsToID115CatalogsAndPlotToPdf(list.of.ID.vcfs, ref.genome = "hg19", 
                                  file = file.path(tempdir(), "VCFToPlotCatID115.counts.test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  catalog.counts.signature <-
    apply(catalog.counts, MARGIN = 2, function(x) x / sum(x))
  catalog.counts.signature <-
    as.catalog.for.ID115(catalog.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <-
    PlotID115CatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatID115.counts.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  if (Sys.getenv("ICAMS.SAVE.TEST.PDF") != "") {
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatID115.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-ID115",
                                 "PlotCatID.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatID115.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-ID115",
                                 "PlotCatID115.counts.signature.test.pdf"))
  } else {
    unlink(file.path(tempdir(), "PlotCatID115.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatID115.counts.signature.test.pdf"))
  }
  graphics.off()
  
})
