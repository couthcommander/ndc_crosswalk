library(readxl)

stdz_ndc <- function(x) {
  sub('^[0]+', '', gsub(' ', '', gsub('-', '', x)))
}

stdz_drug <- function(x) {
  tolower(x)
}

load_ndf <- function(ndc_file) {
  z <- readxl::read_excel(ndc_file, col_types = 'text')
  class(z) <- 'data.frame'
  z1 <- z[,c('NATIONAL_FORMULARY_NAME', 'GENERIC_NAME', 'VA_PRODUCT_NAME', 'VA_CLASSIFICATION_CODE', 'VA_CLASSIFICATION_DESCRIPTION', 'NDC_NUMBER')]
  z1[,'ndc'] <- stdz_ndc(z1[,'NDC_NUMBER'])
  z1[,'drug'] <- stdz_drug(z1[,'VA_PRODUCT_NAME'])
  dat <- unique(z1[order(z1[,'ndc']),c('ndc','drug')])
  row.names(dat) <- NULL
  dat
}
