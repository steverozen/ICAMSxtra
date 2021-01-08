context("AnnotateIDVCFsWithTransRanges")

test_that("AnnotateIDVCFsWithTransRanges function with hg19", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  load("testdata/test_AnnotateIDVCFsWithTransRanges.Rdata")
  list.of.ID.vcfs <- ReadStrelkaIDVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.vcf")
  list <- 
    AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs,
                                  ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
                                  trans.ranges = ICAMS::trans.ranges.GRCh37,
                                  vcf.names = "test.vcf")
  list1 <- AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs, ref.genome = "GRCh37", ICAMS::trans.ranges.GRCh37, vcf.names = "test.vcf")
  list2 <- AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs, ref.genome = "hg19", ICAMS::trans.ranges.GRCh37, vcf.names = "test.vcf")
  expect_equal(list[[1]], annotated.strelka.ID.vcf.GRCh37)
  expect_equal(list, list1)
  expect_equal(list, list2)
})

test_that("AnnotateIDVCFsWithTransRanges with hg38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  load("testdata/test_AnnotateIDVCFsWithTransRanges.Rdata")
  list.of.ID.vcfs <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh38.vcf")
  list3 <- AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs,
                                         ref.genome =
                                           BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                         trans.ranges = ICAMS::trans.ranges.GRCh38,
                                         vcf.names = "test.vcf")
  list4 <- AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs, ref.genome = "GRCh38", ICAMS::trans.ranges.GRCh38, vcf.names = "test.vcf")
  list5 <- AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs, ref.genome = "hg38", ICAMS::trans.ranges.GRCh38, vcf.names = "test.vcf")
  expect_equal(list3[[1]], annotated.strelka.ID.vcf.GRCh38)
  expect_equal(list3, list4)
  expect_equal(list3, list5)
})

