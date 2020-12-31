# Analysis of Screens for Spear-ATAC
# 12/30/20
# High-throughput single-cell chromatin accessibility CRISPR screens enable 
# unbiased identification of regulatory networks in cancer (2020)
# Created by Jeffrey Granja, Sarah Pierce
#
# This script is for analyzing data from Spear-ATAC MS.
#
# NOTE : This script was adapted for better use outside of the original paper. 
# Please report reasonable issues you encounter to the Github and we will 
# try to address them in a reasonable time manner.

###########################################
# Analysis of Spear-ATAC-Screens
###########################################
library(ArchR)
addArchRGenome("hg38")
addArchRThreads(20)

#Helpful functions designed for SpearATAC Analysis
source("Scripts/SpearATAC-Functions.R")

#Input Fragments
#Note we had 6 replicates for K562 but we will show 2 here for simplicity
inputFiles <- c(
	"K562_R1"="data/K562_R1.fragments.tsv.gz", 
	"K562_R2"="data/K562_R2.fragments.tsv.gz"
)

#Get Valid Barcodes
validBC <- getValidBarcodes(
	csvFiles = c("data/K562_R1.singlecell.csv", "data/K562_R2.singlecell.csv"),
	sampleNames = c("K562_R1", "K562_R2")
)

#Get sgRNA Assignments Files that names match the sample names in the ArchRProject
#Files were created from JJJ.R
sgRNAFiles <- c("K562_R1"="data/K562_R1.sgRNA.rds", "K562_R2"="data/K562_R2.sgRNA.rds")

#Create Arrow Files
ArrowFiles <- createArrowFiles(inputFiles=inputFiles, validBarcodes = validBC)
ArrowFiles
# [1] "K562_R2.arrow" "K562_R1.arrow"

#Make an ArchR Project
proj <- ArchRProject(ArrowFiles, outputDirectory = "K562_LS")

#Plot QC of scATAC-seq data
p1 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4)
p2 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "PromoterRatio", plotAs = "violin", alpha = 0.4)
p3 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin", alpha = 0.4)
p4 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin", alpha = 0.4)
plotPDF(p1, p2, p3, p4, name = "QC-singles-scATAC", addDOC = FALSE, width = 4, height = 4)

#Plot Pseudo-Bulk Average QC of scATAC-seq data
p1 <- plotTSSEnrichment(proj)
p2 <- plotFragmentSizes(proj)
plotPDF(p1, p2, name = "QC-pseudobulk-scATAC", addDOC = FALSE, width = 4, height = 4)

#Create an sgRNA assignment matrix
#Iterate over each sgRNA aligned file
sgAssign <- lapply(seq_along(sgRNAFiles), function(x){
	
	message(x)

	#Read in sgRNA Data Frame
	sgDF <- readRDS(sgRNAFiles[x])

	#Create sparseMatrix of sgRNA assignments
	sgMat <- createSpMat(sgDF[,1], sgDF[,2])

	#Create Column Names that match EXACTLY with those in the ArchR Project
	#Cell barcodes in our case were the reverse complement to thos in the scATAC-seq data
	colnames(sgMat) <- paste0(names(sgRNAFiles)[x],"#", reverseComplement(DNAStringSet(colnames(sgMat))),"-1")

	#Check This
	if(sum(colnames(sgMat) %in% proj$cellNames) < 2){
		stop("x=",x,"; Error matching of sgRNA cell barcodes and scATAC-seq cell barcodes was unsuccessful. Please check your input!")
	}

	#Filter those that are in the ArchR Project
	sgMat <- sgMat[,colnames(sgMat) %in% proj$cellNames,drop=FALSE]

	#Compute sgRNA Staistics
	df <- DataFrame(
		cell = colnames(sgMat), #Cell ID
		sgAssign = rownames(sgMat)[apply(sgMat, 2, which.max)], #Maximum sgRNA counts Assignment
		sgCounts =  apply(sgMat, 2, max), #Number of sgRNA counts for max Assignment
		sgTotal = colSums(sgMat), #Number of total sgRNA counts across all Assignments
		sgSpec = apply(sgMat, 2, max) / colSums(sgMat) #Specificity of sgRNA assignment
	)

	#Return this dataframe
	df

}) %>% Reduce("rbind", .)

#Make the rownames the cell barcodes
rownames(sgAssign) <- sgAssign[,1]

#Print
sgAssign
# DataFrame with 10387 rows and 5 columns
#                                              cell    sgAssign  sgCounts
#                                       <character> <character> <numeric>
# K562_R1#TGAATCGTCGGTCCGA-1 K562_R1#TGAATCGTCGGT..  sgBCLAF1-3     30644
# K562_R1#CCCTGATAGCAACTGG-1 K562_R1#CCCTGATAGCAA..     sgTBP-1      6608
# K562_R1#AACTTGGTCTGGTACA-1 K562_R1#AACTTGGTCTGG..    sgNFE2-3      4301
# K562_R1#GATTAGCTCCTCCTGA-1 K562_R1#GATTAGCTCCTC..    sgATF1-3      8893
# K562_R1#AGCTATGGTCCAAGTT-1 K562_R1#AGCTATGGTCCA..    sgPBX2-1      2257
# ...                                           ...         ...       ...
# K562_R2#AATGGCTAGGCAAGGG-1 K562_R2#AATGGCTAGGCA..    sgKLF1-1         2
# K562_R2#GCTCAGGCAATGCACT-1 K562_R2#GCTCAGGCAATG..    sgATF3-1         1
# K562_R2#ATTTGTCAGTACGCGA-1 K562_R2#ATTTGTCAGTAC..    sgZZZ3-1         1
# K562_R2#TCAAGCATCTTAGTGG-1 K562_R2#TCAAGCATCTTA..    sgBRF2-2         1
# K562_R2#ATAGTCGAGCTACGTT-1 K562_R2#ATAGTCGAGCTA..   sgsgNT-12         1
#                              sgTotal    sgSpec
#                            <numeric> <numeric>
# K562_R1#TGAATCGTCGGTCCGA-1     30744  0.996747
# K562_R1#CCCTGATAGCAACTGG-1      6667  0.991150
# K562_R1#AACTTGGTCTGGTACA-1     10500  0.409619
# K562_R1#GATTAGCTCCTCCTGA-1      8915  0.997532
# K562_R1#AGCTATGGTCCAAGTT-1      2271  0.993835
# ...                              ...       ...
# K562_R2#AATGGCTAGGCAAGGG-1         2         1
# K562_R2#GCTCAGGCAATGCACT-1         1         1
# K562_R2#ATTTGTCAGTACGCGA-1         1         1
# K562_R2#TCAAGCATCTTAGTGG-1         1         1
# K562_R2#ATAGTCGAGCTACGTT-1         1         1


#Plot Cutoffs of sgRNA Assignment (NOTE: You may want to adjust these cutoffs based on your results!)
nSg <- 20
Spec <- 0.8
p <- ggPoint(log10(sgAssign$sgCounts), sgAssign$sgSpec, colorDensity = TRUE) +
  geom_vline(xintercept = log10(nSg), lty = "dashed") + 
  geom_hline(yintercept = Spec, lty = "dashed") +
  xlab("log10(sgCounts)") + ylab("Specificity")
plotPDF(p, name = "Plot-Assignment-Density", addDOC = FALSE, width = 5, height = 5)

#Add all this information to your ArchRProject
proj <- addCellColData(proj, data = sgAssign$sgAssign, name = "sgAssign", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgCounts, name = "sgCounts", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgTotal, name = "sgTotal", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgSpec, name = "sgSpec", cells = rownames(sgAssign), force = TRUE)

#Identify those sgRNA Assignemnets that passed your cutoffs
proj$sgAssign2 <- NA
proj$sgAssign2[which(proj$sgCounts > nSg & proj$sgSpec > Spec)] <- proj$sgAssign[which(proj$sgCounts > nSg & proj$sgSpec > Spec)]
proj$sgAssign3 <- stringr::str_split(proj$sgAssign2,pattern="\\-",simplify=TRUE)[,1]

#Print Numbers
table(proj$sgAssign3)
# sgARID2  sgARID3A    sgATF1    sgATF3  sgBCLAF1    sgBRF2     sgCAD   sgCDC5L 
#     122        80       129        64       115       156        69        90 
# sgCEBPB   sgCEBPZ    sgCTCF    sgCUX1    sgELF1   sgFOSL1   sgGABPA   sgGATA1 
#     115       131        99       137       121       147       106       105 
# sgGTF2B   sgHINFP   sgHSPA5    sgKLF1   sgKLF16     sgMAX     sgMYC    sgNFE2 
#      64       123        96        81        83        84        55       129 
#  sgNFYB    sgNRF1    sgPBX2  sgPOLR1D    sgREST    sgRPL9  sgSETDB1    sgsgNT 
#     117       113       148       118       125        92       137       642 
#   sgTBP   sgTFDP1   sgTHAP1  sgTRIM28     sgYY1  sgZBTB11 sgZNF280A  sgZNF407 
#     112       123        93       132       102       116       151       132 
#  sgZZZ3 
#     127 

#########################################################################
#We suggest saving your progress at this moment
#########################################################################
saveRDS(proj, "Save-ArchRProject-W-sgRNA-Assignments.rds")

#We now want to clean our sgAssignments based on the homogeneity of the sgRNA signal. This analysis shouldnt filter more than 5-10%
#of sgRNA assignments. This will help resolve your differential analyses but is not crucial for downstream analysis.
proj <- cleanSgRNA(ArchRProj = proj, groupSg = "sgAssign3", individualSg = "sgAssign2", nonTarget = "sgsgNT")
# ArchR logging to : ArchRLogs/ArchR-cleanSgRNA-17c3e274757af-Date-2020-12-30_Time-21-02-58.log
# If there is an issue, please report to github with logFile!
# Cleaning sgRNA : sgCUX1 (1 of 40)
# TP rate = 0.941605839416058
# Cleaning sgRNA : sgATF3 (2 of 40)
# TP rate = 0.953125
# Cleaning sgRNA : sgBCLAF1 (3 of 40)
# TP rate = 0.991304347826087
# Cleaning sgRNA : sgRPL9 (4 of 40)
# TP rate = 1
# Cleaning sgRNA : sgBRF2 (5 of 40)
# TP rate = 0.967948717948718
# ...

#Clean Up sgRNA assignments based on Purity Ratio (NOTE you may need to adjust this cutoff)
proj$sgAssignFinal <- "UNK"
proj$sgAssignFinal[proj$PurityRatio >= 0.9] <- proj$sgAssign3[proj$PurityRatio >= 0.9]
proj$sgIndividual <- ifelse(proj$sgAssignFinal=="UNK", "UNK", proj$sgAssign2)

#Lets create an unbiased LSI-UMAP to sgRNA assignments by creating an iterativeLSI reduction + UMAP
proj <- addIterativeLSI(proj)
proj <- addUMAP(proj)

#Create Color Palettes
pal1 <- paletteDiscrete(paste0(unique(proj$sgAssign2)))
pal1["NA"] <- "lightgrey"

pal2 <- paletteDiscrete(paste0(unique(proj$sgAssign3)))
pal2["NA"] <- "lightgrey"

pal3 <- paletteDiscrete(paste0(unique(proj$sgAssignFinal)))
pal3["UNK"] <- "lightgrey"

#Plot UMAP Embeddings
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample")
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssign2", pal = pal1)
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssign3", pal = pal2)
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssignFinal", pal = pal3)
plotPDF(p1, p2, p3, p4, name = "Plot-UMAP", addDOC=FALSE)

#Call Peaks using sgAssignments
proj <- addGroupCoverages(proj, groupBy = "sgAssignFinal", force = TRUE)
proj <- addReproduciblePeakSet(proj, 
	groupBy = "sgAssignFinal", force = TRUE, 
	pathToMacs2 = "/Users/jeffreygranja/Library/Python/2.7/bin/macs2")

#Print
getPeakSet(proj)
# GRanges object with 150260 ranges and 12 metadata columns:
#           seqnames              ranges strand |     score
#              <Rle>           <IRanges>  <Rle> | <numeric>
#       UNK     chr1       827241-827741      * |  63.95609
#       UNK     chr1       842709-843209      * |  43.06617
#    sgsgNT     chr1       869692-870192      * |  22.07370
#       UNK     chr1       876691-877191      * |   3.81893
#    sgsgNT     chr1       898606-899106      * |   3.95586
#       ...      ...                 ...    ... .       ...
#    sgNRF1     chrX 155889260-155889760      * |   5.04490
#       UNK     chrX 155893506-155894006      * |  11.81289
#   sgFOSL1     chrX 155944829-155945329      * |   2.71296
#   sgGTF2B     chrX 155966826-155967326      * |   8.36496
#     sgYY1     chrX 155968125-155968625      * |  24.86905

#Create PeakMatrix from Peak Set
proj <- addPeakMatrix(proj, force = TRUE)

#Create Motif Deviations Matrix using Vierstra Motifs
motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
proj <- addMotifAnnotations(proj, motifPWMs = motifPWMs, name = "Vierstra")
proj <- addDeviationsMatrix(proj, peakAnnotation = "Vierstra", force = TRUE)

#Filter These sgRNA non-targetting because they seem to exhibit some differences from the other sgNT
#This step helps a bit, but is not necessary to get differential results
bgd <- grep("sgNT", proj$sgIndividual, value=TRUE) %>% unique
bgd <- bgd[!grepl("-11|-12|-8|-5|-6", bgd)]
proj$sgAssignClean <- proj$sgAssignFinal
proj$sgAssignClean[grepl("-11|-12|-8|-5|-6", proj$sgIndividual)] <- "UNK"

####################################
# Rank Differential Motifs
####################################

#Get Motif Matrix
seMotif <- getMatrixFromProject(proj, "VierstraMatrix")

#Filter Cells
seMotif <- seMotif[,colData(seMotif)$PurityRatio >= 0.9]
seMotif <- seMotif[,!grepl("-11|-12|-8|-5|-6", colData(seMotif)$sgIndividual)] 

#Compute Average Motif Accessibiltiy per Target
gM <- ArchR:::.groupMeans(assays(seMotif)[[2]], colData(seMotif)$sgAssignClean, sparse = TRUE)

#Subtract sgNT
gM_D <- gM %>% {. - .[,grep("sgNT",colnames(.))]}

#Identify Most variable Motif per motif clusters
rMax <- rowMaxs(abs(gM_D))
motifs <- rownames(gM_D)[order(abs(rMax), decreasing = TRUE)] %>% 
	{.[!duplicated(stringr::str_split(.,pattern=":",simplify=TRUE)[,3])]}

#Create Data Frame for Differntials Ranked by Difference
df_D <- reshape2::melt(abs(gM_D[motifs,])) %>% {.[order(.$value,decreasing = TRUE), ]}
df_D$rank <- seq_len(nrow(df_D))
df_D$sgRNA <- stringr::str_split(df_D$Var2,pattern=":",simplify = TRUE)[,1]
df_D <- df_D[df_D$sgRNA %ni% "sgsgNT", ]

#Function for deduplicating after nth appearance
#This is used for labeling to nth motifs per sgRNA
dedupN <- function(x, n){
	idx <- match(x, unique(x))
	idxN <- c()
	tabC <- rep(0, length(unique(x)))
	for(i in seq_along(x)){
		tabC[idx[i]] <- tabC[idx[i]] + 1
		if(tabC[idx[i]] <= n){
			idxN <- c(idxN, i)
		}
	}
	idxN
}

#Identify Top 4 motifs per sgRNA above 0.8
top_D <- df_D[dedupN(df_D$sgRNA, 4),] %>% {.[.$value>0.8,]}

#Create Nice Label
top_D$label <- paste0(top_D$sgRNA,":",stringr::str_split(top_D$Var1,pattern="_",simplify=TRUE)[,1])

#Palette
pal <- paletteDiscrete(df_D$sgRNA)

#Plot
p <- ggplot(df_D[order(df_D$value), ], aes(rank, value, color = sgRNA)) +
	ArchR:::.geom_point_rast2() +
	theme_ArchR() +
	scale_color_manual(values=pal) +
	ylab("Deviation Difference") +
	ggrepel::geom_text_repel(data=top_G, aes(rank, value, label = label))

plotPDF(p, name = "Plot-Top-Motif-Hits-Per-sgRNA", width=6, height=6)

####################################
# Compute Differential Peaks
####################################

#Sort sgRNA Targets so Results are in alphabetical order
useGroups <- sort(unique(proj$sgAssignClean)[unique(proj$sgAssignClean) %ni% c("UNK")])
bgdGroups <- unique(grep("sgNT", proj$sgAssignClean,value=TRUE))

#Differential Peaks
diffPeaks <- getMarkerFeatures(
	ArchRProj = proj, 
	testMethod = "binomial",
	binarize = TRUE,
	useMatrix = "PeakMatrix",
	useGroups = useGroups, 
	bgdGroups = bgdGroups, 
	groupBy = "sgAssignClean",
	bufferRatio = 0.95,
	maxCells = 250,
	threads = 10
)

#Heatmap
heatmapPeaks <- plotMarkerHeatmap(diffPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0")

#Plot Heatmap
plotPDF(heatmapPeaks, name = "Peaks-Marker-Heatmap", width = 8, height = 12, ArchRProj = proj)

#Save RData
save.image("Save-SpearATAC-Analysis.Rdata")

#Session Info
writeLines(sessionInfo(), file = "R-Session-Info.txt")


