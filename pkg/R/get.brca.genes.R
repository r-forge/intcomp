
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti, leo.lahti@iki.fi. All rights reserved.

# Known breast cancer gene list reference
# http://www.nature.com/onc/journal/v18/n56/abs/1203335a.html
# Corresponding file names: tgdb_by_name.cgi.html and tgdb.txt
# downloaded 5.6.2010 from
# http://www.tumor-gene.org/cgi-bin/TGDB/tgdb_by_name.cgi
# and stored to the tgdb object

get.brca.genes <- function (data.genes, xx, tgdb) {

  # Format the symbols
   gen <- toupper(unique(as.character(tgdb[,"Gene"])))
   syn <- toupper(unique(gsub(" ","",unlist(lapply(tgdb[, "Synonyms"], function (x) {strsplit(as.character(x), "\\,")})))))
  symb <- unique(c(gen, syn))
  symb <- unlist(lapply(symb, function (x){strsplit(x, "\\,")[[1]][[1]]}))

  # get Entrez GeneIDs for the gene symbols
  tgdb.genes <- unique(sym2gid(symb, xx))

  # Pick just those cancer genes for which data is available
  cancerGenes <- unique(unlist(intersect(data.genes, tgdb.genes)))

  cancerGenes

}

