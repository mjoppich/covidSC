library(Seurat)
library(stringr)
library(data.table)
library(tidyverse)

library(ggplot2)
library(cowplot)
library("xlsx")

library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/mnt/f/dev/data/Rprojs/covidSC")

readRDS("obj.monocytes_fine_12.Rds")

deResTT_mc_fine = makeDEResults(obj.monocytes_fine_12, assay="RNA", test="t")
exprdfTT_mc_fine = getDEXpressionDF(obj.monocytes_fine_12, deResTT_mc_fine, assay="RNA")
write.table(exprdfTT_mc_fine, "monocytes_fine.de.t.tsv", sep="\t", row.names=F, quote = F)


system("rm analyseMarkers.py")

system("wget https://raw.githubusercontent.com/mjoppich/scrnaseq_celltype_prediction/master/analyseMarkers.py")
system("/usr/bin/python3 analyseMarkers.py --expr-mean mean.cluster --expressing-cell-count anum.cluster --cluster-cell-count num.cluster --organs \"Immune system\" --markers monocytes_fine.de.t.tsv --seurat")

DefaultAssay(obj.monocytes_fine_12) = "RNA"
FeaturePlot(obj.monocytes_fine_12, features = c("FCGR3A", "CD14"), pt.size = 0.01, label=T)
RidgePlot(obj.monocytes_fine_12, features = c("FCGR3A", "CD14"))

#The classical monocyte is characterized by high level expression of the CD14 cell surface receptor
#(CD14++ FCGR3A− monocyte)
# all, except 9,11,3,13
#The non-classical monocyte shows low level expression of CD14 and additional co-expression of the FCGR3A receptor
#(CD14+FCGR3A++ monocyte).
# monocytes_3



FeaturePlot(obj.monocytes_fine_12, features = c("FCGR3A", "CD14"), slot = "data", pt.size = 0.01, label=T)

idents.ncmc = c("monocytes_3")
cells_ncmc  = makeGrpCells(obj.monocytes_fine_12, idents.ncmc)

idents.cmc = setdiff(setdiff(unique(obj.monocytes_fine$sub_cluster_monocytes), idents.ncmc), c("monocytes_11", "monocytes_12", "monocytes_9", "monocytes_3"))
cells_cmc  = makeGrpCells(obj.monocytes_fine_12, idents.cmc)

de.ncmc_cmc = compareCells(obj.monocytes_fine_12, cells_ncmc, cells_cmc, "ncmc", "cmc")

#logfc > 0 => more in ncmc than in cmc

deResList = list("ncmc_cmc"=de.ncmc_cmc)

write.xlsx(de.ncmc_cmc, file = "todoplots/ncmc/de.ncmc_cmc.xlsx" )

saveRDS(deResList, "ncmc_cmc_monocytes.Rds")

#Rscript raAnalysis.R "ncmc_cmc_monocytes.Rds" "go_universe.Rds" "human"
#Rscript goAnalysis.R "ncmc_cmc_monocytes.Rds" "go_universe.Rds"
#Rscript gseAnalysis.R "ncmc_cmc_monocytes.Rds"

ncmc_cmc_monocytes.ra = readRDS("ra.ncmc_cmc_monocytes.Rds")
ncmc_cmc_monocytes.go = readRDS("go.ncmc_cmc_monocytes.Rds")
ncmc_cmc_monocytes.gse = readRDS("gse.ncmc_cmc_monocytes.Rds")


t=makePlotsListsGO(ncmc_cmc_monocytes.go, deResList, "todoplots/ncmc/go/")
t=makePlotsListsGO(ncmc_cmc_monocytes.ra, deResList, "todoplots/ncmc/ra/", name.ggo=NULL, name.ego="rao", name.ego_up="rao_up", name.ego_down="rao_down", name.title="ReactomePA")
t=makePlotsListsGSE(ncmc_cmc_monocytes.gse, deResList, "todoplots/ncmc/gse/")
