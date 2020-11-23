makesummary = function(a, suffix)
{
  out = {}
  out["num"] = length(a)
  
  if (length(a) == 0)
  {
    f = c(0,0,0,0,0)
    meanA = 0
  } else {
    f = fivenum(a)
    meanA = mean(a)
  }
  
  out["min"] = f[1]
  out["lower_hinge"] = f[2]
  out["median"] = f[3]
  out["upper_hinge"] = f[4]
  out["max"] = f[5]
  out["mean"] = meanA
  
  names(out) = paste(names(out), suffix, sep=".")
  
  return(out)
}

getExprData = function(markerObj, markerCells, sampleSuffix, assay="RNA")
{
  expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = "data", assay=assay)
  allgenes = rownames(expTable)
  cellnames = colnames(expTable)
  
  expt.r = as(expTable, "dgTMatrix")
  expt.df = data.frame(r = expt.r@i + 1, c = expt.r@j + 1, x = expt.r@x)
  
  DT <- data.table(expt.df)
  res = DT[, as.list(makesummary(x, sampleSuffix)), by = r]
  res[[paste("anum", sampleSuffix, sep=".")]] = length(cellnames)
  res$gene = allgenes[res$r]
  
  res = res[,r:=NULL]
  
  return(res)
}

getDEXpressionDF = function ( scdata, markers, assay="SCT" )
{
  
  outDF = NULL
  DefaultAssay(object=scdata) = assay  
  clusterIDs = as.character(sort(unique(Idents(scdata))))
  
  scCells = Idents(scdata)
  scCells = names(scCells)
  scCells = unlist(as.character(scCells))
  
  for (clusterID in clusterIDs){
    
    print(clusterID)
    
    cellIdents = Idents(scdata)
    cellIdents.c = names(cellIdents[cellIdents == clusterID])
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))  
    
    cellIdents.bg = setdiff(unlist(lapply(names(cellIdents), as.character)), cellIdents.c)
    
    expvals = getExprData(scdata, cellIdents.c, "cluster", assay=assay)
    expvals.bg = getExprData(scdata, cellIdents.bg, "bg", assay=assay)
    
    modmarkers = markers[[clusterID]]
    modmarkers$gene = rownames(modmarkers)
    
    markerdf = as.data.frame(modmarkers)
    
    if ((nrow(markerdf) > 0) && (nrow(expvals) > 0))
    {
      expvals = merge(markerdf, expvals, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    if ((nrow(expvals) > 0) && (nrow(expvals.bg) > 0))
    {
      expvals = merge(expvals, expvals.bg, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    expvals = as.data.frame(cbind(clusterID, expvals))
    
    if (!is.data.frame(outDF) || nrow(outDF)==0)
    {
      outDF = expvals
    } else {
      outDF = as.data.frame(rbind(outDF, expvals))
    }
    
  }
  
  return(outDF)
  
}

makeDEResults = function(inobj, assay="SCT", test="wilcox")
{
  clusterIDs = as.character(sort(unique(Idents(inobj))))
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
    
    
    cellIdents = Idents(inobj)
    cellIdents.c = names(cellIdents[cellIdents == clusterID])
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))
    
    print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
    
    deMarkers = FindMarkers(inobj, assay=assay, ident.1 = cellIdents.c, test.use=test)
    
    
    retList[[clusterID]] = deMarkers
    
  }
  
  return(retList)
  
}

makeGrpCells = function(scdata, clusterIdents)
{
  cellIdents = Idents(scdata)
  cellIdents.1 = names(cellIdents[cellIdents %in% clusterIdents])
  cellIdents.1 = unlist(lapply(cellIdents.1, as.character))   
  
  return(cellIdents.1)
}


compareCells = function(scdata, cellsID1, cellsID2, suffix1, suffix2, prefix="cluster", test="t", assay="RNA", outfolder="./", all=FALSE)
{
  logfc.threshold = 0.25
  
  if (all==TRUE)
  {
    logfc.threshold = 0.01  
  }
  
  markers = FindMarkers(scdata, assay=assay, ident.1 = cellsID1, ident.2 = cellsID2, test.use=test, logfc.threshold=logfc.threshold)
  
  outvalues1 = getExprData(scdata, cellsID1, suffix1, assay=assay)
  outvalues2 = getExprData(scdata, cellsID2, suffix2, assay=assay) 
  
  
  markers$gene = rownames(markers)
  joinedData = merge(markers, outvalues1, by="gene", all=T)
  joinedData = merge(joinedData, outvalues2, by="gene", all=T)  
  
  joinedData = joinedData[!is.na(joinedData$p_val),]
  
  outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".tsv", sep="")
  
  message(outfile)
  write.table(joinedData, file=outfile, row.names = F,  quote=FALSE, sep='\t')
  
  return(joinedData)
}


makeGSEAnalysis = function(degenes, reverseLogFC=FALSE)
{
  print("Selecting Genes")
  
  # select significantly regulated genes
  siggenes = degenes[degenes$p_val_adj < 0.05,]
  
  if (! "gene" %in% colnames(siggenes))
  {
    siggenes$gene = rownames(siggenes)
  }
  
  # reverse avg logFC
  if (reverseLogFC)
  {
    siggenes$avg_logFC = -siggenes$avg_logFC
  }
  
  
  #siggenes <- siggenes[order(-siggenes$avg_logFC),]
  
  # convert gene symbols to entrez ID for gseGO
  geneNames = bitr(siggenes$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  siggenes = siggenes[siggenes$gene %in% geneNames$SYMBOL,]
  
  geneVector = as.vector(siggenes$avg_logFC)
  names(geneVector) = geneNames$ENTREZID
  # sort decreasing ...
  # we sort by decreasing logFC: geneList = sort(geneList, decreasing = TRUE) https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
  geneVector = sort(geneVector, decreasing = TRUE)
  
  geneVector <- geneVector[!is.na(geneVector)]
  
  print("Got geneVector")
  print(length(geneVector))
  
  if (length(geneVector) == 0)
  {
    return(list(valid=FALSE, gseall=NULL, gsesym=NULL, gseent=NULL))
  }
  
  print("Doing gseGO")
  #enrichGORes <- enrichGO(names(geneVector), ont="BP", OrgDb=org.Hs.eg.db)
  
  # perform gene set enrichment
  gsecc <- gseGO(geneList=geneVector, ont="BP", OrgDb=org.Hs.eg.db, verbose=T, by="fgsea", minGSSize=10, maxGSSize = min(50, length(geneVector)-1), pvalueCutoff=0.5) #nPerm=1000,
  
  print("gseGO Done")
  
  # create a result version sorted by abs(NES) and set size < 50 and qvalue < 0.05
  gseccf = gsecc
  gseccf@result = gseccf@result[order(abs(gseccf@result$NES), decreasing = TRUE),]
  gseccf@result = gseccf@result[gseccf@result$qvalues<0.05 & gseccf@result$setSize<50,]
  
  # create a result version sorted by abs(NES) and set size < 50 and qvalue < 0.05 with GENE SYMBOLS!
  gseccx <- setReadable(gsecc, 'org.Hs.eg.db', 'ENTREZID')
  gseccx@result = gseccx@result[order(abs(gseccx@result$NES), decreasing = TRUE),]
  gseccx@result = gseccx@result[gseccx@result$qvalues<0.05,]
  gseccx@result = gseccx@result[gseccx@result$qvalues<0.05 & gseccx@result$setSize<50,]
  
  return(list(valid=TRUE, gseall=gsecc, gsesym=gseccx, gseent=gseccf))
}

plotGSE = function(clusterID, expression, gseRes, invLFC)
{
  
  cid = as.character(clusterID)
  
  gsePlot = gseRes[[cid]]$gsesym
  
  geneVec = expression[[cid]]$avg_logFC
  
  if (invLFC)
  {
    geneVec = -geneVec
  }
  
  if ("gene" %in% colnames(expression[[cid]]))
  {
    names(geneVec) = expression[[cid]]$gene  
  } else {
    names(geneVec) = rownames(expression[[cid]])
  }
  
  
  p=cnetplot(gsePlot, categorySize="pvalue", foldChange=geneVec, colorEdge=T, showCategory=5) + ggtitle(paste("Cluster", cid, "GSE"))
  
  return(p)
  
}

makeAllGSE = function( deList, reverseLogFC=TRUE )
{
  
  allDEResultIDs = names(deList)
  gseResults = list()
  
  for (resID in allDEResultIDs)
  {
    print(resID)
    
    resDF = deList[[resID]]
    
    print(dim(resDF))
    gseRes = makeGSEAnalysis(resDF, reverseLogFC=reverseLogFC)
    
    gseResults[[resID]] = gseRes
    
  }
  
  return(gseResults)
}