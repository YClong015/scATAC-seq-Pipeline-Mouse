### Function for perforing differential accessibility testing
### groups cells to create pseduo-bulk replicates.
### DA testing performed on pseudo-bulk replicates using DESeq2
apply_DESeq2_test_seurat <- function(seurat.object, 
                                     population.1 = NULL, 
                                     population.2 = NULL, 
                                     exp.thresh = 0.075,
                                     fc.thresh=0.1, 
                                     adj.pval.thresh = 0.05, 
                                     num.splits = 6, 
                                     seed.use = 1234,
                                     replicates.1 = NULL,
                                     replicates.2 = NULL,
                                     verbose = TRUE, 
                                     return.deseq2.res = FALSE, 
                                     assay.use = "ATAC",
                                     ncores = 1) {
  
  if (!'DESeq2' %in% rownames(x = installed.packages())) {
    stop("Please install DESeq2 before using this function
         (http://bioconductor.org/packages/release/bioc/html/DESeq2.html)")
  }
  
  if (is.null(population.1) & is.null(replicates.1)) {
    stop("Both population.1 and replicates.1 cannot be NULL")
  }
  
  if (!is.null(population.1) & !is.null(replicates.1)) {
    stop("Values for either population.1 OR replicates.1 should be provided - not both")
  }
  
  ## reduce counts in a cluster to num.splits cells
  if (is.null(replicates.1)) {
    high.expressed.peaks <- GetExpressedPeaks(seurat.object, population.1, population.2, threshold = exp.thresh, assay.use = assay.use)
    length(high.expressed.peaks)
  } else {
    high.expressed.peaks <- GetExpressedPeaks(seurat.object, unlist(replicates.1), unlist(replicates.2), threshold = exp.thresh, assay.use = assay.use)
  }

  if (verbose) print(paste(length(high.expressed.peaks), "accesible peaks passing threshold ", toString(exp.thresh)))
  
  peaks.use <- high.expressed.peaks
  if (verbose) print(paste(length(peaks.use), "individual peak sites to test"))
  
  ## make pseudo-bulk profiles out of cells
  ## set a seed to allow replication of results
  set.seed(seed.use)
  if (is.null(replicates.1)) {
    
    if (length(population.1) == 1) {
      cells.1 <- names(Seurat::Idents(seurat.object))[which(Seurat::Idents(seurat.object) == population.1)]
    } else{
      cells.1 <- population.1
    }
    
    cells.1 = sample(cells.1)
    cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets1 <- replicates.1
  }
  
  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- Seurat::GetAssayData(seurat.object, layer = "counts", assay = assay.use)[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- peaks.use
  colnames(profile.set1) <- paste0("Population1_", 1:length(cell.sets1))
  
  ## create a profile set for second cluster
  if (is.null(replicates.2)) {
    if (is.null(population.2)) {
      cells.2 <- setdiff(colnames(seurat.object), cells.1)
    } else {
      if (length(population.2) == 1) {
        cells.2 <- names(Seurat::Idents(seurat.object))[which(Seurat::Idents(seurat.object) == population.2)]
      } else {
        cells.2 <- population.2
      }
    }
    
    cells.2 = sample(cells.2)
    cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets2 <- replicates.2
  }
  
  
  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- Seurat::GetAssayData(seurat.object, layer = "counts", assay = assay.use)[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set2[, i] <- this.profile
    } else {
      profile.set2[, i] <- sub.matrix
    }
  }
  rownames(profile.set2) <- peaks.use
  colnames(profile.set2) <- paste0("Population2_", 1:length(cell.sets2))
  
  ## merge the count matrices together
  peak.matrix <- cbind(profile.set1, profile.set2)
  
  ## Create the sample table
  sampleTable <- data.frame(row.names = c(colnames(profile.set1), colnames(profile.set2)),
                            condition = c(rep("target", ncol(profile.set1)),
                                          rep("comparison", ncol(profile.set2))))
  
  ## Build the DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=peak.matrix, 
                                        colData=sampleTable, 
                                        design=~condition)
  
  ## Run DESeq2
  if (verbose) print("Running DESeq2 test...")
  
  dds <- DESeq2::DESeq(dds)
  res.specific <- DESeq2::results(dds, contrast=c( "condition", "target", "comparison" ), test="Wald")
  
  if (return.deseq2.res) {
    return(dds)
  } else {
    return(res.specific)
  } 

}



GetExpressedPeaks <- function(seurat.object, population.1, population.2=NULL, threshold=0.05, assay.use = "ATAC") {
  
  if (length(population.1) == 1){ # cluster identity used as input
    foreground.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)==population.1])
  } else { # cell identity used as input
    foreground.set = population.1
  }
  if (is.null(population.2)) {
    remainder.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)!=population.1])
  } else {
    if (length(population.2) == 1) { # cluster identity used as input
      remainder.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)==population.2])
    } else { # cell identity used as input
      remainder.set = population.2
    }
  }
  
  peak.names = rownames(seurat.object)
  
  # Get the peaks expressed in the foreground set based on proportion of non-zeros
  this.data <- Seurat::GetAssayData(seurat.object, layer = "counts", assay=assay.use)
  nz.row.foreground = tabulate(this.data[, foreground.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  peaks.foreground = peak.names[which(nz.prop.foreground > threshold)]
  
  # Now identify the peaks expressed in the background set
  nz.row.background = tabulate(this.data[, remainder.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.background = nz.row.background/length(remainder.set)
  peaks.background = peak.names[which(nz.prop.background > threshold)]
  
  return(union(peaks.foreground, peaks.background))
}
