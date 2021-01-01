suppressPackageStartupMessages(library(Biostrings)) #Bioconductor
suppressPackageStartupMessages(library(readr)) #CRAN
suppressPackageStartupMessages(library(data.table)) #CRAN

nucleotideSubstitution <- function(str, letters = c("A","C","G","T","N"), n = 1){
  temp <- paste0(str)
  strx <- unlist(lapply(seq_len(nchar(temp[[1]])), function(y){
    message(y)
    unlist(lapply(seq_along(letters), function(z){
      tempz <- temp
      substring(tempz, y, y) <- letters[z]
      tempz
    }))
  }))
  dt <- data.table(new = strx, original = rep(temp, nchar(temp[[1]]) * length(letters)))
  dt <- dt[order(dt[,2]),]
  dt <- DataFrame(data.frame(dt))
  dt[,2] <- Rle(dt[,2])
  gc()
  dt
}

#1bp mismatch
barcodes <- DNAStringSet(readr::read_tsv("10x_Barcodes.txt", col_names = FALSE)$X1)
barcodes <- reverseComplement(barcodes)
barcodesMisMatch <- nucleotideSubstitution(paste0(barcodes), n = 1)
saveRDS(barcodesMisMatch, "10x-scATAC-Barcodes-1-MisMatch-190613.rds")
