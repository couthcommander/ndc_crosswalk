##################################################
# package requirements
##################################################

library(DBI)
library(RSQLite) # if using default `dbCon`
library(data.table)

##################################################
# configuration
##################################################

# folder with RXNORM downloads -- should contain zip(s), ie `RxNorm_full_12092005.zip`
dataPath <- file.path('zips')
# should RXNORM data be loaded? set to TRUE if new data
reloadData <- FALSE
# function to define database connection -- default is recommended
dbCon <- function() DBI::dbConnect(RSQLite::SQLite(), dbname = 'rxnorm.sqlite')
# name of database table -- default is recommended
dataTable <- 'rxnorm_raw'
# function to write file of bad output; NDC codes with discrepancies
badOutput <- function(data) write.csv(data, file = file.path('ndc_discrep.csv'), row.names = FALSE)
# function to write file of good output; NDC crosswalk
goodOutput <- function(data) write.csv(data, file = file.path('ndc_gs_xwalk.csv'), row.names = FALSE)
# drugname description lookup -- code will likely fail if changed
lu <- read.csv(file.path('rx_desc_lookup.csv'))

##################################################
# utility functions
##################################################

stdz_ndc <- function(x) {
  sub('^[0]+', '', gsub(' ', '', gsub('-', '', x)))
}

stdz_drug <- function(x) {
  tolower(x)
}

bigWord <- function(x) {
  vapply(strsplit(x, "\\s|,|/|\\(|\\)|-"), function(i) {
    i[order(nchar(i), decreasing = TRUE)[1]]
  }, character(1))
}

##################################################
# database functions
##################################################

writeTable <- function(db, table, dat, ...) {
  DBI::dbWriteTable(conn = db, name = table, value = dat, row.names = FALSE, ...)
}

deleteSource <- function(db, table, src) {
  ts <- DBI::dbListTables(db)
  if(!(table %in% ts)) {
    stop('table not found')
  }
  if(grepl('[^-0-9]', src)) {
    stop('src should be YYYY-MM-DD format')
  }
  qry <- sprintf("DELETE FROM %s WHERE src = '%s'", table, src)
  res <- DBI::dbSendQuery(db, qry)
  DBI::dbClearResult(res)
}

tblCheck <- function(db, table) {
  ts <- DBI::dbListTables(db)
  if(!(table %in% ts)) {
    stop('table not found')
  }
  qry <- sprintf("SELECT src, count(*) AS cnt FROM %s GROUP BY src", table)
  res <- DBI::dbSendQuery(db, qry)
  out <- DBI::dbFetch(res)
  DBI::dbClearResult(res)
  out
}

datSource <- function(db, table, src) {
  ts <- DBI::dbListTables(db)
  if(!(table %in% ts)) {
    stop('table not found')
  }
  if(grepl('[^-0-9]', src)) {
    stop('src should be YYYY-MM-DD format')
  }
  qry <- sprintf("SELECT * FROM %s WHERE src = '%s'", table, src)
  res <- DBI::dbSendQuery(db, qry)
  out <- DBI::dbFetch(res)
  DBI::dbClearResult(res)
  out
}

datNDC <- function(db, table, code) {
  ts <- DBI::dbListTables(db)
  if(!(table %in% ts)) {
    stop('table not found')
  }
  isV <- length(code) > 1
  if(isV) {
    if(any(grepl('[^0-9a-zA-Z]', code))) {
      stop('code should be alphanumeric format')
    }
    code <- paste(code, collapse = "', '")
    qry <- sprintf("SELECT * FROM %s WHERE ndc IN ('%s') ORDER BY ndc, src", table, code)
  } else {
    if(grepl('[^0-9a-zA-Z]', code)) {
      stop('code should be alphanumeric format')
    }
    qry <- sprintf("SELECT * FROM %s WHERE ndc = '%s' ORDER BY src", table, code)
  }
  res <- DBI::dbSendQuery(db, qry)
  out <- DBI::dbFetch(res)
  DBI::dbClearResult(res)
  out
}

datQry <- function(db, qry) {
  res <- DBI::dbSendQuery(db, qry)
  out <- DBI::dbFetch(res)
  DBI::dbClearResult(res)
  out
}

##################################################
# RXNORM extraction functions
##################################################

# `match` returns first match
# name often has multiple options matching rxcui

load_rxnorm <- function(sat_file, conso_file, version = 1) {
  sat <- fread(sat_file, sep='|', quote = '')
  sat <- sat[, .(V1, V4, V9, V10, V11, V12)]
  setnames(sat, c('RXCUI','RXAUI','ATN','SAB','ATV','SUPPRESS'))
  sat <- sat[SAB == 'VANDF']
  rxcui <- unique(sat[['RXAUI']])
  if(version == 1) {
    nameOpt <- sat[ATN == 'NF_NAME']
    gnrcOpt <- sat[ATN == 'VA_GENERIC_NAME']
    clssOpt <- sat[ATN == 'VA_CLASS_NAME']
    ndcCode <- unique(sat[ATN == 'NDC', .(RXAUI, ATV,SUPPRESS)])
  } else if(version == 2) {
    nameOpt <- sat[ATN == 'TRN']
    gnrcOpt <- sat[ATN == 'PRN']
    clssOpt <- sat[ATN == 'VMO'] #or VAC?
  }
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
  dat <- unique(dat[order(dat[,'ndc']),c('ndc','drug','SUPPRESS')])
# Suppressible flag. Values = N, O, Y, or E.
# N - not suppressible.
# O - Specific individual names (atoms) set as Obsolete because the name is no longer provided by the original source.
# Y - Suppressed by RxNorm editor.
# E - unquantified, non-prescribable drug with related quantified, prescribable drugs.
  row.names(dat) <- NULL
  dat
}

get <- function(zipfile) {
  td <- tempdir()
  unzip(zipfile, exdir = td)
  if('rrf' %in% list.files(td)) {
    wd <- file.path(td, 'rrf')
  } else {
    wd <- td
  }
  dat <- load_rxnorm(file.path(wd, 'RXNSAT.RRF'), file.path(wd, 'RXNCONSO.RRF'))
  unlink(td, recursive = TRUE)
  dat
}

save2db <- function(path, rawTable, dbc) {
  zips <- list.files(path, pattern='zip$')
  fd <- as.Date(gsub('[^0-9]', '', zips), format = '%m%d%Y')
  names(zips) <- fd
  dats <- zips[order(fd)]
  for(i in seq_along(dats)) {
    dat <- get(file.path(path, dats[[i]]))
    src <- names(dats)[i]
    dat <- cbind(dat, src = src)
    deleteSource(dbc, rawTable, src)
    writeTable(dbc, rawTable, dat, append = TRUE)
  }
  src
}

##################################################
# data de-duplication functions
##################################################

find_dup <- function(x) {
  i <- x[,'ndc']
  x[i %in% i[duplicated(i)],]
}

stdzdrug <- function(x) {
  x <- sub('_#[0-9]', '', tolower(x))
  x[x == 'a & d oint'] <- 'vitamin a/vitamin d oint'
  x <- sub('apap ', 'acetaminophen ', x)
  x <- sub('asa ', 'aspirin ', x)
  x <- sub('hctz ', 'hydrochlorothiazide', x)
  x <- sub('\\(alc-f)[ ]?', '', x)
  x <- sub('\\((p|s)f)[ ]?', '', x)
  x <- sub('heparin na \\(pork)', 'heparin na', x)
  x <- sub('globulin,human', 'globulin', x)
  x <- sub('cream,top', 'cream', x)
  x <- sub('tab[, ]ud', 'tab', x)
  x <- sub('tab,sa,ud', 'tab,sa', x)
  x <- sub('inj kit', 'inj,kit', x)
  x <- sub('inj[, ]soln', 'inj', x)
  x <- sub('(inj,conc|conc inj)', 'inj', x)
  x <- sub('mg/vil', 'mg/ml', x)
  x <- sub('unt/', 'unit/', x)
  x <- sub('na succ ', 'na succinate ', x)
  x <- sub('triamcinolone (acetonide|acet)', 'triamcinolone', x)
  x <- sub('(midazolam|dexmedetomidine) hcl ', '\\1 ', x)
  x <- sub('(naloxone) hcl[ ]*', '\\1 ', x)
  x <- sub('(naloxone) \\(eqv-narcan)[ ]*', '\\1 ', x)
  x <- sub('soln,itrc', 'inj', x)
  # ignore all numbers
  x <- gsub('[ ]*[0-9][0-9,]*[ ]*', ' X', x)
  x <- sub('inj, Xml amp', 'inj', x)
  x <- sub('inj,amp, Xml', 'inj', x)
#   x <- sub('([^/]{10,})/.*', '\\1', x)
  x
}

snp <- function(x) {
  paste(sort(x), collapse = ' ')
}

diffs <- function(x, y) {
  xd <- stdzdrug(x[,'drug'])
  yd <- stdzdrug(y[,'drug'])
  xd <- vapply(strsplit(xd, '[ /,;]+'), snp, character(1))
  yd <- vapply(strsplit(yd, '[ /,;]+'), snp, character(1))
  ix <- match(x[,'ndc'], y[,'ndc'])
  sd <- stringdist::stringdist(xd, yd[ix], method = 'jw')
  fix <- which(sd > 0)
  fsd <- sd[fix]
  data.frame(ndc = x[fix,'ndc'], drug.x = xd[fix], drug.y = yd[ix[fix]], fsd)[order(fsd),]
}

checkdup <- function(dat) {
  dat1 <- find_dup(dat)
  dat1 <- dat1[order(dat1[,'ndc']),]
  x <- dat1[,'ndc']
  nx <- length(x)
  dx <- duplicated(x)
  x_s <- which(!dx)
  key <- seq(nx) - rep(x_s, diff(c(x_s, nx + 1))) + 1
  d3 <- diffs(dat1[key==1,], dat1[key==2,])
  d3x <- bigWord(d3[,'drug.x'])
  d3y <- bigWord(d3[,'drug.y'])
  d3[,'maineq'] <- +(d3x == d3y)
  d3
}

diffByLabel <- function(x, labels) {
  # `x` should have "d1" and "d2"
  a <- x
  a[,c('label1','label2')] <- NA_character_
  nrowA <- nrow(a)
  nrowL <- nrow(labels)
  # `labels` should have "rxname" and "desc"
  commonCombos <- c('moisturizer, skin irritation', 'decongestant, pain relief', 'allergy, skin irritation',
  'nutritional supplement, minerals', 'nutritional supplement, eye vitamins', 'nutritional supplement, moisturizer',
  'nutritional supplement, antioxidant', 'nutritional supplement, probiotic',
  'ace inhibitor, diuretic', 'laxative, minerals', 'laxative, incontinence', 'dry eye, probiotic',
  'vaccine-hib, vaccine-hepb, vaccine-meningitis', 'skin irritation, hemorrhoid',
  'antacid, aspirin', 'antacid, gas relief', 'barrier, moisturizer', 'barrier, antiseptic',
  'syringe, incontinence', 'syringe, insulin'
  )
  ud <- sort(c(unique(labels[,'desc']), commonCombos))
  combos <- grep(',', ud)
  cmbgrp <- strsplit(ud[combos], ', ')
  lc <- length(combos)
  for(i in seq(lc)) {
    cmbgrp[[i]] <- c(cmbgrp[[i]], ud[combos[i]])
  }
  inCombo <- function(x) {
    nx <- length(x)
    tmp <- vapply(x, function(j) vapply(cmbgrp, function(i) j %in% i, logical(1)), logical(lc))
    methA <- combos[rowSums(tmp) == nx]
    methB <- combos[ud[combos] %in% x]
    union(methA, methB)
  }
  inCombo1 <- function(x) {
    inCombo(x)[1]
  }
  dco <- match(labels[,'desc'], ud)
  m <- matrix(0, nrowA*2, length(ud))
  src <- sprintf("%s([^a-zA-Z]|$)", labels[,'rxname'])
  for(i in seq(nrowL)) {
    i1 <- grep(src[i], a[,'d1'])
    i2 <- grep(src[i], a[,'d2'])
    m[i1, dco[i]] <- 1
    m[nrowA + i2, dco[i]] <- 1
  }
  ix <- which(rowSums(m) > 1)
  for(i in seq_along(ix)) {
    spot <- which(m[ix[i],] > 0)
    opt <- inCombo1(ud[spot])
    if(!is.na(opt) && length(opt) > 0) {
      m[ix[i],] <- 0
      m[ix[i], opt] <- 1
    }
  }

  for(i in seq(nrowA)) {
    a1 <- ud[m[i,] == 1]
    a2 <- ud[m[nrowA + i,] == 1]
    if(length(a1)) a[i,'label1'] <- a1
    if(length(a2)) a[i,'label2'] <- a2
  }

  # allow combos to match
  for(i in seq(nrowA)) {
    nl1 <- a[i,'label1']
    nl2 <- a[i,'label2']
    if(nl1 == nl2) next
    nlc1 <- inCombo(nl1)
    nlc2 <- inCombo(nl2)
    if(!is.na(nlc1[1])) {
      nl1 <- ud[nlc1]
    }
    if(!is.na(nlc2[1])) {
      nl2 <- ud[nlc2]
    }
    inAll <- intersect(nl1, nl2)
    if(length(inAll) == 0) {
      inAll <- intersect(unlist(strsplit(nl1, ', ')), unlist(strsplit(nl2, ', ')))
      if(length(inAll) > 1) {
        inAll <- paste(inAll, collapse = ', ')
      }
    }
    if(length(inAll) == 1) {
      a[i,'label1'] <- inAll
      a[i,'label2'] <- inAll
    }
  }
  a[,'labeleq'] <- +(a[,'label1'] == a[,'label2'])
  a
}

ndc_compete <- function(codes, tbl, dbc) {
  info <- datNDC(dbc, tbl, codes)
  rxl <- tapply(info[,'drug'], info[,'ndc'], unique)
  rxm <- mapply(expand.grid, rxl, rxl, SIMPLIFY = FALSE, stringsAsFactors = FALSE)
  for(i in seq_along(rxm)) {
    rxm[[i]] <- cbind(names(rxm)[i], t(apply(rxm[[i]], 1, sort)))
  }
  rxd <- do.call(rbind, c(rxm, make.row.names = FALSE))
  df <- as.data.frame(unique(rxd[rxd[,2] != rxd[,3],]))
  names(df) <- c('ndc','d1','d2')
  rownames(df) <- NULL
  df
}

##################################################
# create output
##################################################

con <- dbCon()
if(reloadData) {
  save2db(dataPath, dataTable, con)
}

qry <- sprintf('select ndc, drug from %s order by ndc, src', dataTable)
allndc <- datQry(con, qry)
dat <- allndc[rev(!duplicated(rev(allndc[,'ndc']))),]
ub <- unique(allndc)
ubc <- checkdup(ub)
tubc <- table(ub[,'ndc'])
dupndcs <- names(tubc[tubc > 2])
dupndc <- c(setdiff(ubc[ubc[,'maineq'] == 0,'ndc'], dupndcs), dupndcs)
dd <- ndc_compete(dupndc, dataTable, con)

isDup <- dat[,'ndc'] %in% dd[,'ndc']
xw1 <- dat[!isDup,]

zz <- diffByLabel(dd, lu)
ndc_ne <- zz[zz[,'labeleq'] == 0, 'ndc']
ndc_eq <- zz[zz[,'labeleq'] == 1, 'ndc']
if(length(ndc_ne)) {
  d1 <- datNDC(con, dataTable, ndc_ne)
  dup_bad <- d1[!duplicated(do.call(paste, c(d1[,c('ndc','drug')], sep = '|'))), c('ndc','drug','src')]
  badOutput(dup_bad)
}
if(length(ndc_eq)) {
  d2 <- datNDC(con, dataTable, ndc_eq)
  xw2 <- d2[rev(!duplicated(rev(d2[,'ndc']))), c('ndc','drug')]
  xw <- rbind(xw1, xw2)
  xw <- xw[order(xw[,'ndc']),]
  goodOutput(xw)
} else {
  goodOutput(xw1)
}

DBI::dbDisconnect(con)
