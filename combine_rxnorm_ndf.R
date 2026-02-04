source('load_ndf.R')

# source('load_rxnorm.R')
# rx <- load_rxnorm('RXNSAT.RRF', 'RXNCONSO.RRF')
## rather than load from single RXNORM data set, use `rxnorm_ndc.R` script
rx <- read.csv('ndc_gs_xwalk.csv')
ndf <- load_ndf('PharmacyProductSystem_NationalDrugCodeExtract.csv', skip = 2)

# 10-digit vs 11-digit
# 00000-0000-00
# 00000-0000-0
# 00000-000-00
# 0000-0000-00

findName <- function(x, v) {
  x[grep(v, x[,'drug']),]
}

findNDC <- function(x, v) {
  x[grep(v, x[,'ndc']),]
}

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

diffs <- function(x, y) {
  xd <- stdzdrug(x[,'drug'])
  yd <- stdzdrug(y[,'drug'])
  ix <- match(x[,'ndc'], y[,'ndc'])
  sd <- stringdist::stringdist(xd, yd[ix], method = 'jw')
  fix <- which(sd > 0)
  fsd <- sd[fix]
  list(
    x[is.na(ix),],
    cbind(x[fix,], drug.y = y[ix[fix],'drug'], fsd)[order(fsd),]
  )
}

snp <- function(x) {
  paste(sort(x), collapse = ' ')
}

diffs2 <- function(x, y) {
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

## step 1: deduplicate data set

## no more duplicates
find_dup(rx)
find_dup(ndf)

## step 2: compare discrepancies

com_ids <- intersect(ndf[,'ndc'], rx[,'ndc'])
d1 <- rx[match(com_ids, rx[,'ndc']),]
d2 <- ndf[match(com_ids, ndf[,'ndc']),]
view_diff1 <- diffs(d1, d2)
view_diff2 <- diffs2(d1, d2)

## step 3: add together

add_ndf <- ndf[match(setdiff(ndf[,'ndc'], rx[,'ndc']), ndf[,'ndc']),]
mylist <- rbind(rx, add_ndf)
mylist <- mylist[order(mylist[,'ndc']),]
write.csv(mylist, file = 'ndc_xwalk.csv', row.names = FALSE)
