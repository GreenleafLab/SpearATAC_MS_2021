message("Loading Pacakges...")

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("ShortRead", quietly = TRUE)) BiocManager::install("ShortRead")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("fastmatch", quietly = TRUE)) install.packages("fastmatch")

suppressPackageStartupMessages(library(ShortRead)) #Bioconductor
suppressPackageStartupMessages(library(optparse)) #CRAN
suppressPackageStartupMessages(library(readr)) #CRAN
suppressPackageStartupMessages(library(magrittr)) #CRAN
suppressPackageStartupMessages(library(data.table)) #CRAN
suppressPackageStartupMessages(library(fastmatch)) #CRAN
suppressPackageStartupMessages(library(Biostrings)) #Bioconductor

######################################################################
######################################################################
######################################################################

nucleotideSubstitution <- function(str, letters = c("A","C","G","T","N"), n = 1){
  temp <- paste0(str)
  n <- nchar(temp[[1]])
  strx <- unlist(lapply(seq_len(nchar(temp[[1]])), function(y){
    message(y, " of ", n, " Combinations Generated...")
    unlist(lapply(seq_along(letters), function(z){
      tempz <- temp
      substring(tempz, y, y) <- letters[z]
      tempz
    }))
  }))
  message("Finished creating barcodes mismatches...")
  dt <- data.table(new = strx, original = rep(temp, nchar(temp[[1]]) * length(letters)))
  dt <- dt[order(dt[,2]),]
  dt <- DataFrame(data.frame(dt))
  dt[,2] <- Rle(dt[,2])
  gc()
  dt
}

chunkDictMatch <- function(x, dict, chunkSize = 100000){
  chunks <- split(seq_along(x), ceiling(seq_along(x)/chunkSize))
  matchIdx <- lapply(seq_along(chunks), function(i){
    idx <- which(vcountPDict(dict, x[chunks[[i]]]) > 0, arr.ind = TRUE)
    idx[,2] <- chunks[[i]][idx[,2]]
    idx
  }) %>% Reduce("rbind", .)
  matchIdx
}

######################################################################
######################################################################
######################################################################
option_list <- list(
  make_option(c("--sgRNA"), type="character", default="sgRNA-LargeScreen.txt"),
  make_option(c("--barcodes"), type="character", default="10x_Barcodes.txt"),
  make_option(c("--fastq_read"), type="character", default="fastq-Read1"),
  make_option(c("--fastq_index"), type="character", default="fastq-Index2"),
  make_option(c("--output"), type="character", default = "Save-Alignment.rds")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
######################################################################
######################################################################
######################################################################

#guideRNA
message("Reading in sgRNA...")
sgDF <- data.frame(readr::read_tsv(opt$sgRNA, col_names = FALSE))
sgPDict_forward <- PDict(paste0(DNAStringSet(sgDF[,2])), tb.start = 1, tb.end = min(width(BStringSet(sgDF[,2]))))
sgPDict_reverse <- PDict(paste0(reverseComplement(DNAStringSet(sgDF[,2]))), tb.start = 1, tb.end = min(width(BStringSet(sgDF[,2]))))

message("Creating mismatches...")
barcodes <- DNAStringSet(readr::read_tsv(opt$barcodes, col_names = FALSE)$X1)
barcodes <- reverseComplement(barcodes)
barcodesMisMatch <- nucleotideSubstitution(paste0(barcodes), n = 1)

#Create Read Streams
message("Creating read streams...")
chunkSize <- 5*10^5
readStream   <- FastqStreamer(opt$fastq_read,  n = chunkSize)
indexStream  <- FastqStreamer(opt$fastq_index, n = chunkSize)

alignedSG <- DataFrame()
start <- Sys.time()
n <- 0
message("Demultiplexing and adding Barcodes...")
repeat{

  message("Processed ", round(n/10^6,6), " M fastq reads, Found ", 
    round(nrow(alignedSG)/10^6,6), " M reads, ", round(difftime(Sys.time(),start,units="mins"),2), " minutes...")
  
  #Yield Fastq Info
  fqread <- yield(readStream)
  fqindex <- yield(indexStream)
  
  if(length(fqread) == 0 | length(fqindex) == 0){
    break
  }

  #Get Read Names
  readname <- BStringSet(id(fqread))
  index2name <- BStringSet(id(fqindex))

  if(identical(readname,index2name)){
    
    #ID ok sgRNA
    afreq <- matrixStats::rowMaxs(alphabetFrequency(sread(fqread))) < min(width(sread(fqread)))
    
    if(length(afreq) > 0){
      
      matchIdx <- fmatch(paste0(sread(fqindex[afreq])), barcodesMisMatch[,1])
      idx2 <- !is.na(matchIdx)

      if(length(which(idx2)) > 0){

        #Time to Align
        sgID <- sread(fqread)[afreq][idx2]
        cellID <- BStringSet(barcodesMisMatch[matchIdx[idx2],2])

        #Forward
        sgAlign <- DataFrame(chunkDictMatch(sgID, sgPDict_forward)) #Make Sure Memory Doesnt Go Nuts
        if(nrow(sgAlign) > 0){
          sgAlign[,2] <- cellID[sgAlign[,2]]
          sgAlign[,3] <- "Forward"
          alignedSG <- rbind(alignedSG, sgAlign)
        }

        #Reverse
        sgAlign <- DataFrame(chunkDictMatch(sgID, sgPDict_reverse)) #Make Sure Memory Doesnt Go Nuts
        if(nrow(sgAlign) > 0){
          sgAlign[,2] <- cellID[sgAlign[,2]]
          sgAlign[,3] <- "Reverse"
          alignedSG <- rbind(alignedSG, sgAlign)
        }

      }

    }

  }else{
    stop("Error readname ids do not match!!!!")
  }

  n <- n + length(fqread)

}

#Convert to sgRNA names
alignedSG[,1] <- sgDF[alignedSG[,1],1]

#Save Output
saveRDS(alignedSG, opt$output)





