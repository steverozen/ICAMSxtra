context("Categorizing the mutation type of deletions (into 115 categories)")

test_that("CreateOneColID115Matrx insertions", {
  MakeTestInsVCF <- function() {
    return(data.frame(
      seq.context = c("TTTTTTTTTTTTCGACCCCCCCCCCCC",
                      "TTTTTTTTTTTTGAACCCCCCCCCC",
                      "TTTTTTTTTTTTGCCCCCCCCCCCC",
                      "TTTTTTTTTTTTGCCCCCCCCCCCC"),
      REF = c("C", "G", "G", "G"),
      ALT = c("CGA", "GA", "GA", "GC"),
      seq.context.width = c(12, 12, 12, 12),
      stringsAsFactors = FALSE
    ))
  }
  load("testdata/create_one_col_insert_test_ID115.Rdata")
  expect_equal(
    CreateOneColID115Matrix(MakeTestInsVCF(), NULL)[[1]],
    create.one.col.insert.test.ID115)
})

test_that("CreateOneColID115Matrix deletions", {
  MakeTestDelVCF <- function() {
    return(
      data.frame(
        seq.context = c(
          "GAGGTATACATTGTGTTTACTTTTTCTATGTTTATGTACAATAGTAATATCTTTATAGTTATACTAACGTTATTAAAATAAGTAATTATATTAACTAAGTTTAGGACCAGTTTCTAGT",
          "GACCACTGAGAACCCAGGTTTTAGGCCCACCCCGGTACCAGGCCAGCCCCTGT",
          "AAGGTTTGGCTTCA",
          "ATTAAAATGGGGTT"),
        REF = c("ATAGTTATAC", "GCCCA", "TG", "AT"),
        ALT = c("A", "G", "T", "A"),
        seq.context.width = c(54, 24, 6, 6),
        stringsAsFactors = FALSE))
  }
  load("testdata/create_one_col_delete_test_ID115.Rdata")
  expect_equal(
    CreateOneColID115Matrix(MakeTestDelVCF(), NULL)[[1]],
    create.one.col.delete.test.ID115)
})

