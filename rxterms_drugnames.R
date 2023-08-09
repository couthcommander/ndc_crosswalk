library(data.table)

snp <- function(x) {
  paste(sort(x), collapse = '^')
}

tDrug <- function(name) {
  vapply(strsplit(name, '[/]'), snp, character(1))
}

rxterm <- function(path) {
  zips <- list.files(path, pattern='zip$')
  fd <- gsub('[^0-9]', '', zips)
  names(zips) <- fd
  dats <- zips[order(fd)]
  zipfile <- file.path(path, dats[[length(dats)]])

  td <- tempdir()
  unzip(zipfile, exdir = td)
  dat_file <- list.files(td, pattern = 'RxTerms[0-9]+.txt', full.names = TRUE)
  dat <- data.table::fread(dat_file, sep='|', quote = '')
  unlink(td, recursive = TRUE)

  a <- dat[['DISPLAY_NAME']]
  b <- dat[['ROUTE']]
  pat <- sprintf('[ ]+\\((%s)\\)', paste(unique(b), collapse = '|'))
  d <- tolower(gsub(pat, '', a))
  d <- gsub('[0-9.]+ mg\\/[0-9.]+ ml[ ]+', '', d)
  d <- gsub('[0-9./ ]+(day pack|day|carton|sample pack|mixed pack|kit)$', '', d)
  d <- gsub('[0-9.mg ]+\\/[0-9.mgl ]+(dose pack|taper)$', '', d)
  d <- gsub('[0-9.mg ]+\\/[0-9.mgl ]+$', '', d)
  d <- gsub('[ ]+\\((usp|[0-9]+)\\)', '', d)
  d <- gsub('[ ]+\\([0-9]+ (count|unt)\\)', '', d)
  d <- gsub('[ ]+\\([0-9.]+ mg[^)]*\\)', '', d)
  d <- gsub('[ ]+\\([0-9]+ x [^)]*\\)', '', d)
  d <- gsub('[ ]+\\(level [0-9]+\\)', '', d)
  d <- gsub('[ ]+\\(for patients[^)]+\\)', '', d)
  ix <- grepl('/', d) & !grepl('[0-9 ]/[0-9 ]', d)
  d1 <- tDrug(d[ix])
  d2 <- d[!ix]
  d3 <- unlist(strsplit(d[ix], '[/]'))
  sort(union(d3, d2))
}

checkname <- function(x) {
  paste(rxt[vapply(rxt, grepl, logical(1), x)], collapse = '^')
}

rxt <- rxterm('zip_terms')
# ignore short words
rxt <- rxt[nchar(rxt) >= 4]
addl_rxt <- c(
'amprenavir',
'albumin',
'biperiden',
'bitolterol',
'cantharidin',
'cefditoren',
'cyclizine',
'deoxycorticosterone',
'deracoxib',
'detomidine',
'diatrizoate',
'dihydrotachysterol',
'dipivefrin',
'dyphylline',
'ergonovine',
'firocoxib',
'florfenicol',
'flunixin',
'flurandrenolide',
'fomivirsen',
'gonadorelin',
'guanabenz',
'hylan polymers',
'iodixanol',
'iohexol',
'iopanoic acid',
'iopromide',
'ipodate',
'iso-sulfan blue',
'methicillin sodium',
'milbemycin oxime',
'monascus purpureus west',
'moricizine',
'nitrofurazone',
'octoxynol',
'plicamycin',
'ritodrine',
'sermorelin',
'sibutramine',
'sulfanilamide',
'tadenan',
'tacrine',
'thiopental',
'tocainide',
'tolazoline',
'triclocarban',
'urokinase'
)

rxt <- sort(c(rxt, addl_rxt))
writeLines(rxt, 'rxterms_drugnames.txt')
