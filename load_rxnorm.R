library(data.table)

stdz_ndc <- function(x) {
  sub('^[0]+', '', gsub(' ', '', gsub('-', '', x)))
}

stdz_drug <- function(x) {
  tolower(x)
}

# `match` returns first match
# name often has multiple options matching rxcui

load_rxnorm <- function(sat_file, conso_file, suppress = FALSE) {
  sat <- fread(sat_file, sep='|', quote = '')
  sat <- sat[, .(V1, V4, V9, V10, V11, V12)]
  setnames(sat, c('RXCUI','RXAUI','ATN','SAB','ATV','SUPPRESS'))
  if(suppress) {
    sat <- sat[SAB == 'VANDF' & SUPPRESS == 'N']
  } else {
    sat <- sat[SAB == 'VANDF']
  }
  rxcui <- unique(sat[['RXAUI']])
  nameOpt <- sat[ATN == 'NF_NAME']
  gnrcOpt <- sat[ATN == 'VA_GENERIC_NAME']
  clssOpt <- sat[ATN == 'VA_CLASS_NAME']
  ndcCode <- unique(sat[ATN == 'NDC', .(RXAUI, ATV)])
  setnames(ndcCode, 'ATV', 'NDC')
  name <- nameOpt[match(rxcui, nameOpt[['RXAUI']]), ATV]
  gnrc <- gnrcOpt[match(rxcui, gnrcOpt[['RXAUI']]), ATV]
  clss <- clssOpt[match(rxcui, clssOpt[['RXAUI']]), ATV]

  conso <- fread(conso_file, sep='|', quote = '')
  conso <- conso[, .(V1, V8, V12, V15)]
  setnames(conso, c('RXCUI','RXAUI','SAB','STR'))
  strOpt <- conso[SAB == 'VANDF']
#   strOpt <- conso[SAB == 'RXNORM']
  str <- strOpt[match(rxcui, strOpt[['RXAUI']]), STR]

  res <- data.frame(RXAUI = rxcui, NF_NAME = name, VA_GENERIC_NAME = gnrc, STR = str, VA_CLASS_NAME = clss, stringsAsFactors = FALSE)
  dat <- merge(res, ndcCode)
  dat[,'ndc'] <- stdz_ndc(dat[,'NDC'])
  dat[,'drug'] <- stdz_drug(dat[,'STR']) # NF_NAME, VA_GENERIC_NAME, STR
  dat <- unique(dat[order(dat[,'ndc']),c('ndc','drug')])
  row.names(dat) <- NULL
  dat
}
