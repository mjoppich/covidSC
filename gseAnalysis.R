library(clusterProfiler)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)

inFile = args[1]
universeFile = args[2]

inName = tools::file_path_sans_ext(basename(inFile))
baseDir = dirname(inFile)

outFile = paste(baseDir, paste("gse", inName, "Rds", sep="."), sep="/")

print(outFile)


makeGSEAnalysis = function(degenes, universe, reverseLogFC=FALSE)
{
  print("Selecting Genes")
  
  if ((is.null(degenes)) || (nrow(degenes) == 0))
  {
    print("NULL genes")
    
    return(list(valid=FALSE, gseall=NULL, gsesym=NULL, gseent=NULL))
  }
  
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

  # perform gene set enrichment
  gsecc <- gseGO(geneList=geneVector, ont="BP",
                 OrgDb=org.Hs.eg.db, verbose=T,
                 by="DOSE", nPerm=1000,
                 minGSSize=min(5, max(1,length(geneVector)-1)),
                 #maxGSSize = min(50, max(1, length(geneVector)-1)),
                 pvalueCutoff=0.5) #nPerm=1000,
  
  print("gseGO Done")
  
  # create a result version sorted by abs(NES) and set size < 50 and qvalue < 0.05
  gseccf = gsecc
  gseccf@result = gseccf@result[order(abs(gseccf@result$NES), decreasing = TRUE),]
  #gseccf@result = gseccf@result[gseccf@result$qvalues<0.05 & gseccf@result$setSize<50,]
  
  # create a result version sorted by abs(NES) and set size < 50 and qvalue < 0.05 with GENE SYMBOLS!
  gseccx <- setReadable(gsecc, 'org.Hs.eg.db', 'ENTREZID')
  gseccx@result = gseccx@result[order(abs(gseccx@result$NES), decreasing = TRUE),]
  #gseccx@result = gseccx@result[gseccx@result$qvalues<0.05,]
  #gseccx@result = gseccx@result[gseccx@result$qvalues<0.05 & gseccx@result$setSize<50,]
  
  
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

makeAllGSE = function( deList, universe, reverseLogFC=TRUE )
{
  
  allDEResultIDs = names(deList)
  gseResults = list()
  
  for (resID in allDEResultIDs)
  {
    print(resID)
    
    resDF = deList[[resID]]
    
    print(dim(resDF))
    gseRes = makeGSEAnalysis(resDF, universe, reverseLogFC=reverseLogFC)
    
    gseResults[[resID]] = gseRes
    
  }
  
  return(gseResults)
}

print(inFile)
deList = readRDS(inFile)

outList = makeAllGSE(deList, NULL, reverseLogFC = F)

saveRDS(outList, outFile)
