# Source this file from ICAMSxtra top level directory.

cat(getwd(), "\n")

list.of.ID.vcfs1 <- 
  ICAMS::ReadStrelkaIDVCFs("tests/testthat/testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.vcf")
annotated.strelka.ID.vcf.GRCh37 <- 
  AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs1, ref.genome = "GRCh37", 
                                ICAMS::trans.ranges.GRCh37, vcf.names = "test.vcf")[[1]]

list.of.ID.vcfs2 <- 
  ICAMS::ReadStrelkaIDVCFs("tests/testthat/testdata/Strelka.ID.GRCh38.vcf")
annotated.strelka.ID.vcf.GRCh38 <- 
  AnnotateIDVCFsWithTransRanges(list.of.ID.vcfs2, ref.genome = "GRCh38", 
                                ICAMS::trans.ranges.GRCh38, vcf.names = "test.vcf")[[1]]

save(annotated.strelka.ID.vcf.GRCh37, annotated.strelka.ID.vcf.GRCh38,
    file = "tests/testthat/testdata/test_AnnotateIDVCFsWithTransRanges.Rdata")

