library(singleCellTK)
library(celda)
# Import data
sampleNames <- c("BAL_1", "BAL_4", "Blood_1", "Blood_4")
sce.full <- importCellRanger(cellRangerDirs = rep("../", length(sampleNames)), sampleDirs = sampleNames)
rownames(sce.full) <- rowData(sce.full)$feature_name
sce.full <- scater::logNormCounts(sce.full)

# Remove duplicate gene names
counts(sce.full) <- as.matrix(counts(sce.full))
gene.table <- table(rownames(sce.full))
gene.duplicates <- gene.table[gene.table > 1]
gene.duplicates.names <- names(gene.duplicates)
empty <- character(0)
if (!identical(empty, gene.duplicates.names)){
  for (genename in gene.duplicates.names){
    genename <- gsub(" (1 of many)", "", genename, fixed=TRUE)
    indices <- which(grepl(genename, rownames(sce.full)))
    num <- length(indices)
    for (i in 1:num){
      rownames(sce.full)[indices[i]] <- paste0(genename, "-", i)
    }
  }
}

### QC
sce.full <- runCellQC(sce.full, sample = sce.full$sample, mitoPrefix = "^mt-", mitoGeneLocation = "rownames", algorithms = c("QCMetrics", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX", "decontX_bg"))
reportCellQC(sce.full)
plotRunPerCellQCResults(sce.full, sample = sce.full$sample)

### Filtering
cutoff.nUMI <- 500
cutoff.nGenes <- 500
cutoff.pctMito <- 5
sce <- subsetSCECols(sce.full, colData = c(paste0("total > ",cutoff.nUMI), paste0("detected >", cutoff.nGenes), paste0("mito_percent <", cutoff.pctMito)))
plotRunPerCellQCResults(sce, sample = sce$sample)

K <- 15
L <- 150
sce <- reportCeldaCGRun(sce, K=K, L=L, sampleLabel = sce$sample, output_sce_prefix = paste0("celda_cg_K", K, "_L", L), output_file=paste0("CeldaCG_RunReport_K", K, "_L", L))
neutrophils <- c("Itgam", "Ptprc", "Adgre1", "Ly6g", "Itgb2", "Csf1r")
FeaturePlot(data, features = neutrophils, reduction = "celda_UMAP")
plotDimReduceCluster(sce, reducedDimName = "celda_UMAP")

########################### Removing non-neutrophil clusters
sce <- sce[,which(!(celdaClusters(sce) %in% c(1,8,11,12,13,14,15)))]

K <- 8
L <- 100
sce <- reportCeldaCGRun(sce, K=K, L=L, sampleLabel = sce$sample, output_sce_prefix = paste0("celda_cg_K", K, "_L", L), output_file=paste0("CeldaCG_RunReport_K", K, "_L", L))
sce$clusters <- celdaClusters(sce)
reducedDim(sce, "celda_UMAP") <- reducedDim(altExp(sce), "celda_UMAP")
sce <- findMarkerDiffExp(sce, cluster = celdaClusters(sce))

########################### Additional downstream analysis
## Slingshot 
sce.nest <- sce
sce.nest <- scater::runPCA(sce.nest)
reducedDim(sce.nest, "UMAP") <- reducedDim(altExp(sce.nest), "celda_UMAP")
reducedDim(sce.nest, "TSNE") <- reducedDim(altExp(sce.nest), "celda_tSNE")
sce.sling <- slingshot(sce.nest, reducedDim='PCA', clusterLabels = celdaClusters(sce.nest))
pseudo.paths <- slingPseudotime(sce.sling)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
gg <- scater::plotUMAP(sce.sling, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=celda_UMAP1, y=celda_UMAP2), size=1.2)
}
gg

## RNA Velocity
# import loom file
ldat <- ReadVelocity(file = "/rprojectnb/scscore/Analysis/TraberK/2019-03-22/RNA_velocity/combined_samples.loom")
names <- make.unique(rownames(ldat$spliced))
rownames(ldat$spliced) <- names
rownames(ldat$unspliced) <- names
rownames(ldat$ambiguous) <- names
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# Velocity for neutrophils only
cols <- gsub("-","_",substr(colnames(sce),1, nchar(colnames(sce)) - 2))
bm2 <- bm[,which(colnames(bm) %in% cols)]
sce <- sce[,which(cols %in% colnames(bm2))]
rd <- reducedDim(altExp(sce), "celda_UMAP")
rownames(rd) <- gsub("-","_",substr(rownames(rd),1,nchar(rownames(rd))-2))
identical(rownames(rd), colnames(bm2))
bm2[["UMAP"]] <- CreateDimReducObject(embeddings = rd, key = "UMAP_", assay = DefaultAssay(bm2))
Idents(bm2) <- celdaClusters(sce)
bm2$celltype <- celdaClusters(sce)
bm2 <- RunVelocity(object = bm2, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm2)))
names(x = ident.colors) <- levels(x = bm2)
cell.colors <- ident.colors[Idents(object = bm2)]
names(x = cell.colors) <- colnames(x = bm2)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm2, reduction = "UMAP"), vel = Tool(object = bm2, 
                                                                                              slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


# Velocity for BAL neutrophils only
sce <- sce[,which(sce$sample %in% c("BAL_1", "BAL_4"))]
sce$original_clusters <- celdaClusters(sce)
sce <- sce[,which(!(celdaClusters(sce) %in% c(1,2)))]
celdaClusters(sce) <- factor(celdaClusters(sce))
sce$celda_clusters <- celdaClusters(sce)
reducedDim(sce, "UMAP") <- reducedDim(altExp(sce), "celda_UMAP")
reducedDim(sce, "TSNE") <- reducedDim(altExp(sce), "celda_tSNE")
cols <- gsub("-","_",substr(colnames(sce),1, nchar(colnames(sce)) - 2))
bm3 <- bm[,which(colnames(bm) %in% cols)]
sce <- sce[,which(cols %in% colnames(bm3))]
rd <- reducedDim(altExp(sce), "celda_UMAP")
rownames(rd) <- gsub("-","_",substr(rownames(rd),1,nchar(rownames(rd))-2))
identical(rownames(rd), colnames(bm3))
bm3[["UMAP"]] <- CreateDimReducObject(embeddings = rd, key = "UMAP_", assay = DefaultAssay(bm3))
Idents(bm3) <- celdaClusters(sce)
bm3$celltype <- celdaClusters(sce)
bm3 <- RunVelocity(object = bm3, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm3)))
names(x = ident.colors) <- levels(x = bm3)
cell.colors <- ident.colors[Idents(object = bm3)]
names(x = cell.colors) <- colnames(x = bm3)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm3, reduction = "UMAP"), vel = Tool(object = bm3, 
                                                                                              slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


## VAM Enrichment of custom genesets 
Genesets <- read_csv("Custom_Genesets.csv")
genesets <- lapply(Genesets, function(x) x[!is.na(x)])
sce <- runFindMarker(sce, cluster = "celda_clusters", log2fcThreshold = 0.25)
de <- metadata(sce)$findMarker
genesets$Degranulation <- unique(genesets$Degranulation) # Remove duplicates 
genesets$`Response to Cytokines and Chemokines` <- unique(genesets$`Response to Cytokines and Chemokines`)
genesets$`Response to Bacterial Factors` <- unique(genesets$`Response to Bacterial Factors`)
sce <- importGeneSetsFromList(sce, geneSetList = genesets, collectionName = "CustomGeneSets")
sce <- runVAM(sce, geneSetCollectionName = "CustomGeneSets", useAssay = "logcounts")
for (i in 1:length(genesets)){
  print(plotPathway(sce, resultName = "VAM_CustomGeneSets_CDF", geneset = names(genesets)[i], groupBy = "celda_clusters"))
}


## GO Enrichment
databases <- enrichR::listEnrichrDbs()$libraryName
dbs <- c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023") # Select databases
for (i in levels(as.factor(celdaClusters(sce)))){
  genes <- de[which(de$findMarker_cluster == i), ]$Gene
  sce <- runEnrichR(inSCE = sce, features = genes, db = dbs,
                    analysisName = "tutorial_analysis") 
  res <- getEnrichRResult(sce, analysisName = "tutorial_analysis")$result
  res$Term <- substr(res$Term, 1, stringr::str_length(res$Term) - 13) # Remove GO ID from term names
  res$P.value <- res$Adjusted.P.value # Use adjusted p values for plotting 
  print(plotEnrich(res, numChar = 100, y = "Ratio") + ggtitle(paste("Cluster",i,"GO Enrichment")) +theme(text = element_text(size=12)) + guides(fill = guide_colorbar(title = "Adjusted P")))
}

