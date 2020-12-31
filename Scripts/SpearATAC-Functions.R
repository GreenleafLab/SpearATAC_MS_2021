#Access ArchR Hidden Functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}

#Function for creating a sparse confusion matrix
createSpMat <- function(i,j){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#Function For cleaning sgRNA
cleanSgRNA <- function(
  ArchRProj = NULL, 
  groupSg = NULL,
  individualSg = NULL,
  nonTarget = "sgNT", 
  useMatrix = "TileMatrix",
  scaleTo = 10000, 
  excludeChr = "chrM",
  verbose = FALSE,
  k = 20,
  seed = 1,
  tR  = 0.9,
  varFeatures = 10000,
  dimsToUse = 1:30,
  totalFeatures = 500000,
  filterQuantile = 0.999,
  plotDir = "cleanSgRNA",
  threads = getArchRThreads(),
  logFile = createLogFile("cleanSgRNA")
  ){

  dir.create(plotDir)
  time <- make.names(paste0(Sys.time()))

  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
      tryCatch({
          eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
      }, error = function(x) {
      })
  }
  
  tstart <- Sys.time()

  .startLogging(logFile=logFile)

  #1 Get Info for LSI

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]

  #Check if Matrix is supported and check type
  stopifnot(any(tolower(useMatrix) %in% c("tilematrix","peakmatrix")))
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }
  if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }

  chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)
  
  #Compute Row Sums Across All Samples
  .logDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  if(useMatrix == "TileMatrix"){
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
    totalAcc$start <- (totalAcc$idx - 1) * tileSize
  }else{
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
  }
  gc()

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = "nFrags")
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(e){
      tryCatch({
        .getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
      }, error = function(y){
        stop("Could not determine depth from nFrags or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  #Filter Chromosomes
  if(length(excludeChr) > 0){
    totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
  }

  #2 Create sgRNA Group Matrix
  .logDiffTime("Computing Top Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  rmTop <- floor((1-filterQuantile) * totalFeatures)
  topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), totalFeatures + rmTop)[-seq_len(rmTop)]
  topFeatures <- totalAcc[sort(topIdx),]
  
  #Need to determine which cells are being held out
  sgRNA <- getCellColData(ArchRProj, groupSg)
  groupList <- SimpleList(split(rownames(sgRNA), paste0(sgRNA[,1])))

  groupMat <- .getGroupMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures, 
    groupList = groupList, 
    threads = threads, 
    useIndex = FALSE, 
    verbose = verbose, 
    useMatrix = useMatrix, 
    asSparse = FALSE, 
    tstart = tstart
  )
  groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo

  #Get Groups for individual SgRNA
  sgRNA_Individual <- getCellColData(ArchRProj, individualSg)
  groupIndividualList <- SimpleList(split(rownames(sgRNA_Individual), paste0(sgRNA_Individual[,1])))
  
  indMat <- .getGroupMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures, 
    groupList = groupIndividualList, 
    threads = threads, 
    useIndex = FALSE, 
    verbose = verbose, 
    useMatrix = useMatrix, 
    asSparse = FALSE, 
    tstart = tstart
  )
  indMat <- t(t(indMat) / colSums(indMat)) * scaleTo

  #Unique Sg
  uniqSg <- unique(sgRNA[,1])
  uniqSg <- uniqSg[!is.na(uniqSg)]
  uniqSg <- uniqSg[uniqSg %ni% nonTarget]

  mapSg <- split(paste0(sgRNA_Individual[,1]), paste0(sgRNA[,1]))
  mapSg2 <- lapply(seq_along(mapSg), function(x) unique(mapSg[[x]]))
  names(mapSg2) <- names(mapSg)


  #########################################################
  #
  # Individual sgRNA
  #
  #########################################################

  assignList <- lapply(seq_along(uniqSg), function(x){

    message(sprintf("Cleaning sgRNA : %s (%s of %s)", uniqSg[x], x, length(uniqSg)))

    gM <- groupMat[, c(uniqSg[x], nonTarget)]
    gM <- edgeR::cpm(gM, log = TRUE, prior.count = 0.5)

    iM <- indMat[, mapSg2[[uniqSg[x]]],drop=FALSE]
    iM <- edgeR::cpm(iM, log = TRUE, prior.count = 0.5)

    diffUp <- order(gM[,1] - gM[,2], decreasing = TRUE)
    diffUpIdx <- which(rowSums(iM - gM[,2] > 0) > 1)
    featuresUp <- head(diffUp[diffUp %in% diffUpIdx], floor(varFeatures/2))

    diffDo <- order(gM[,2] - gM[,1], decreasing = TRUE)
    diffDoIdx <- which(rowSums(gM[,2] - iM > 0) > 1)
    featuresDo <- head(diffDo[diffDo %in% diffDoIdx], floor(varFeatures/2))

    variableFeatures <- topFeatures[sort(c(featuresUp, featuresDo)),]

  #1. Get Partial Matrix With All Cells of Interest
  mat <- .getPartialMatrix(
    ArrowFiles = ArrowFiles,
    featureDF = variableFeatures,
    useMatrix = useMatrix,
        cellNames = c(groupList[[uniqSg[x]]], groupList[[nonTarget]]),
    doSampleCells = FALSE,
    threads = threads,
    verbose = FALSE
  )  

  #2. Create Initial Manifold
  LSI_T_NT <- .computeLSI(
    mat = mat[,c(groupList[[uniqSg[x]]], groupList[[nonTarget]])], 
    LSIMethod = 2, 
    scaleTo = scaleTo,
    nDimensions = max(dimsToUse),
    binarize = TRUE, 
    verbose = FALSE, 
    tstart = tstart
  )

  #3. Compute UMAP Manifold
  umap_T_NT <- uwot::umap(
    X = .rowZscores(LSI_T_NT[[1]]), 
    n_neighbors = 40, 
    min_dist = 0.4, 
    metric = "cosine", 
    ret_model = TRUE
  )

  #Convert To data.frame
  umap_T_NT <- data.frame(umap_T_NT[[1]])
  rownames(umap_T_NT) <- rownames(LSI_T_NT[[1]])
  colnames(umap_T_NT) <- c("UMAP1", "UMAP2")

  #Compute KNN
  knnIdx <- .computeKNN(umap_T_NT, umap_T_NT, k = k)
  rownames(knnIdx) <- rownames(umap_T_NT)
  targetRatio <- rowSums(knnIdx <= length(groupList[[uniqSg[x]]])) / k

  TPT <- intersect(names(which(targetRatio >= tR)), groupList[[uniqSg[x]]])
  TPNT <- intersect(names(which(targetRatio <= 1 - tR)), groupList[[nonTarget]])

  message(sprintf("TP rate = %s", length(TPT) / length(groupList[[uniqSg[x]]])))

  knnIdx2 <- .computeKNN(umap_T_NT[c(TPT,TPNT),], umap_T_NT, k = k)
  rownames(knnIdx2) <- rownames(umap_T_NT)
  targetRatio2 <- rowSums(knnIdx2 <= length(TPT)) / k

  p0 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], targetRatio, discrete=F,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p1 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], targetRatio2, discrete=F,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p2 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], paste0(seq_len(nrow(umap_T_NT)) <= length(groupList[[uniqSg[x]]])), discrete=T,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
    pdf(file.path(plotDir, paste0(uniqSg[x], ".pdf")), width = 8, height = 4)
    ggAlignPlots(p0, p1, p2, type="h", draw=TRUE)
    dev.off()

    df <- data.frame(cellNames = names(targetRatio2), targetRatio = targetRatio2, sgRNA = sgRNA[names(targetRatio2),])

  SimpleList(targetRatio = df, UMAP = umap_T_NT)

  })

  df2 <- lapply(seq_along(assignList), function(x){
  assignList[[x]][[1]]
  }) %>% Reduce("rbind", .)

  df2 <- df2[order(df2[,2], decreasing=TRUE),]
  df2 <- df2[!duplicated(paste0(df2[,1])), ]

  sgRNA$PurityRatio <- as.numeric(rep(-1, nrow(sgRNA)))
  sgRNA[paste0(df2[,1]),"PurityRatio"] <- df2[,2]

  sgRNA[which(paste0(sgRNA[,1]) == nonTarget),2] <- 1 - sgRNA[which(paste0(sgRNA[,1]) == nonTarget),2]
  colnames(sgRNA) <- c("sgRNA", "PurityRatio")

  ArchRProj <- addCellColData(ArchRProj, data = sgRNA[,"PurityRatio"], name = "PurityRatio", cells = rownames(sgRNA))

  ArchRProj

}

