library(clusterProfiler)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)

inFile = args[1]
universeFile = args[2]

inName = tools::file_path_sans_ext(basename(inFile))
baseDir = dirname(inFile)

outFile = paste(baseDir, paste("go", inName, "Rds", sep="."), sep="/")

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
  
  #geneRegVec = siggenes$avg_logFC
  #names(geneRegVec) = siggenes$gene
  
  
  # convert gene symbols to entrez ID for gseGO
  geneNames = bitr(siggenes$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  print("Unresolved gene symbols")
  print(siggenes[!(siggenes$gene %in% geneNames$SYMBOL),]$gene)
  
  siggenes = siggenes[siggenes$gene %in% geneNames$SYMBOL,]
  
  geneVector = as.vector(siggenes$avg_logFC)
  names(geneVector) = geneNames$ENTREZID
  # sort decreasing ...
  # we sort by decreasing logFC: geneList = sort(geneList, decreasing = TRUE) https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
  geneVector = sort(geneVector, decreasing = TRUE)
  
  #geneVector = names(geneVector)
  geneVector <- geneVector[!is.na(geneVector)]
  
  print("Got geneRegVec")
  print(length(geneVector))
  
  
  if (length(geneVector) == 0)
  {
    return(list(valid=FALSE, ggo=NULL, ego=NULL))
  }
  
  print(paste(length(geneVector), paste(head(geneVector), collapse=";")))
  print(head(universe))
  
  ggo <- groupGO(gene     = names(geneVector),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)
  
  ego = enrichGO(gene=names(geneVector),
                 universe=universe,
                 ont="BP",
                 OrgDb=org.Hs.eg.db,
                 readable=TRUE,
                 minGSSize=5,
                 pvalueCutoff=0.5,
                 qvalueCutoff=0.5,
                 pAdjustMethod = "BH")
  upGenes=geneVector[geneVector > 0]
  print(paste("UpGenes", length(upGenes)))
  
  ego_up = enrichGO(gene=names(upGenes),
                 universe=universe,
                 ont="BP",
                 OrgDb=org.Hs.eg.db,
                 readable=TRUE,
                 minGSSize=5,
                 pvalueCutoff=0.5,
                 pAdjustMethod = "BH")
  
  downGenes = geneVector[geneVector < 0]
  print(paste("downGenes", length(downGenes)))
  
  ego_down = enrichGO(gene=names(downGenes),
                 universe=universe,
                 ont="BP",
                 OrgDb=org.Hs.eg.db,
                 readable=TRUE,
                 minGSSize=5,
                 pvalueCutoff=0.5,
                 pAdjustMethod = "BH")
  
  print(head(ego))
  
  return(list(valid=TRUE, ggo=ggo, ego=ego, ego_up=ego_up, ego_down=ego_down))
}


makeAllGSE = function( deList, universe, reverseLogFC=TRUE )
{
  
  allDEResultIDs = names(deList)
  gseResults = list()
  
  universeSym = universe
  geneNames = bitr(universeSym, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  universeEntrez = geneNames[universeSym %in% geneNames$SYMBOL,]$ENTREZID
  
  universeEntrez = universeEntrez[!is.na(universeEntrez)]
  
  for (resID in allDEResultIDs)
  {
    print(resID)
    
    if (resID != "nk.asympt_sympt")
    {
      #next
    }

    resDF = deList[[resID]]
    
    print(dim(resDF))
    gseRes = makeGSEAnalysis(resDF, universeEntrez, reverseLogFC=reverseLogFC)
    
    gseResults[[resID]] = gseRes
    
  }
  
  return(gseResults)
}

deList = readRDS(inFile)
universe = readRDS(universeFile)

outList = makeAllGSE(deList, universe, reverseLogFC = F)

saveRDS(outList, outFile)
