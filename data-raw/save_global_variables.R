# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/create_catalog_row_headers_ID115.R")

usethis::use_data(catalog.row.headers.ID115,
                  internal = TRUE,
                  overwrite = TRUE)