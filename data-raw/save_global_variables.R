# Source this file from ICAMSxtra top level directory.

cat(getwd(), "\n")

source("data-raw/create_catalog_row_headers_ID115.R")
source("data-raw/create_catalogs_ID115.R")
source("data-raw/create_trans_proportions.R")
source("data-raw/create_testdata_for_AnnotateIDVCFsWithTransRanges.R")
source("data-raw/create_spectra_from_vcfs.R")

usethis::use_data(catalog.row.headers.ID115,
                  internal = TRUE,
                  overwrite = TRUE)

usethis::use_data(catalog.row.order,
                  GRCh37.proportions,
                  GRCh38.proportions,
                  overwrite = TRUE)