context("PlotID115Catalog")

test_that("PlotID115Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(2, 2, 2, 1))
  on.exit(par(opar))
  test.file <- c("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.vcf")
  list.of.ID.vcfs <- ReadStrelkaIDVCFs(test.file)
  list.of.catalog <- VCFsToID115Catalogs(list.of.ID.vcfs, ref.genome = "hg19",
                                 region = "genome")
  catalog <- list.of.catalog$catalog
  out <- PlotID115Catalog(catalog)
  out1 <- PlotID115Catalog(catalog, ylim = c(0, 1000))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  cat.counts.signature <- apply(catalog, MARGIN = 2, function(x) x / sum(x))
  cat.counts.signature <-
    as.catalog.for.ID115(cat.counts.signature, ref.genome = "hg19",
               region = "genome", catalog.type = "counts.signature")
  out <- PlotID115Catalog(cat.counts.signature)
  out1 <- PlotID115Catalog(cat.counts.signature, ylim = c(0, 0.5))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  graphics.off()
})