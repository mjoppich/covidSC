

library(Tempora)
library(Seurat)
library(RCurl)
library(tidyverse)
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(scatterpie)
library(RColorBrewer)
library(igraph)
library(ggrepel)
library(stringr)
library(scales)

save_plot = function(plotobj,outname, fig.width, fig.height)
{
  print(paste(outname, fig.width, fig.height))
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()
  
  return(plotobj)
}

if (!require('RCurl')) {
  install.package('RCurl')
} 


PlotVaryingPWs <- function(object, identifier, outprefix, pw.max=50){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@varying.pws)){
    print("IdentifyVaryingPWs has not been run or no temporally varying pathways were detected. Please run IdentifyVaryingPWs or re-run with a more relaxed p-value cutoff See ?Tempora::IdentifyVaryingPWs for details")
    return()
  }

  varying_pathways <- sort(object@varying.pws)

  if (length(varying_pathways) > pw.max)
  {
    varying_pathways = varying_pathways[1:pw.max]
  }
  gsva_bycluster <- object@cluster.pathways
  gams <- object@gams
    
    annotDF = object@cluster.metadata[, c("Id", "label")]
    annotDF$chrID = as.character(annotDF$Id)
    clusterLabels = annotDF$label
    clusterLabels = unlist(lapply(str_split(annotDF$label, "-"), function(x){x[2]}))

    names(clusterLabels) = annotDF$chrID


  cat("\nPlotting time-dependent pathways...")
  print(paste("Number of time-dependent pathways", length(varying_pathways)))

  if (!dir.exists(outprefix))
  {
    dir.create(outprefix)
  }

  for (i in 1:length(varying_pathways)){
    fname=paste(outprefix, paste("PlotVaryingPWs", identifier, i, sep="_"), sep="/")
    print(paste("Saving to file", fname))
    plot.new()
    png(filename=paste(fname, "png", sep="."), width = 8, height = 6, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)

      par(mar=c(5,5,5, 10))      
        
      if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) > 1){
        plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]))
        plot_df$time <- object@cluster.metadata$Cluster_time_score
        plot_df$cluster = clusterLabels[plot_df$cluster]
      }
      else if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) == 1) {
        plot_df <- data.frame(cluster=names(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ])
        plot_df$time <- object@cluster.metadata$Cluster_time_score
        plot_df$cluster = clusterLabels[plot_df$cluster]
      }
      id <- which(names(gams)==names(varying_pathways)[i])
      mgcv::plot.gam(gams[[id[1]]], main = paste0(names(varying_pathways)[i]), xlab = "Inferred time", ylab="Pathway expression level", bty="l", cex.main = 1, xaxt = "n", shade= F, se=3, scheme=1)
      xmin <- par("usr")[1]
      xmax <- par("usr")[2]
      points(x=plot_df$time, y=plot_df$value, pch=20, col="navy", cex=0.9, xpd=NA)
      text(x=plot_df$time, y=plot_df$value, labels=plot_df[,1], pos = 4, cex = 0.5, col="navy", xpd=NA)
      legend("topright", legend = "Cluster", pch = 20, col = "navy", bty="n", text.col="navy", cex=0.9)
      legend("topright", legend=paste0("\nAdjusted p-value = ", round(varying_pathways[[i]], 5)), bty="n", cex=0.9)
      axis(side=1, at=c(xmin, xmax), labels = c("Early", "Late"), tick=T)

    dev.off()


    
    print(paste("Saving to file", fname))
    pdf(file=paste(fname, "pdf", sep="."),width = 8, height = 6)

      par(mar=c(5,5,5, 10))      
        
      if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) > 1){
        plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]))
        plot_df$time <- object@cluster.metadata$Cluster_time_score
        plot_df$cluster = clusterLabels[plot_df$cluster]
      }
      else if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) == 1) {
        plot_df <- data.frame(cluster=names(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ])
        plot_df$time <- object@cluster.metadata$Cluster_time_score
        plot_df$cluster = clusterLabels[plot_df$cluster]
      }
      id <- which(names(gams)==names(varying_pathways)[i])
      mgcv::plot.gam(gams[[id[1]]], main = paste0(names(varying_pathways)[i]), xlab = "Inferred time", ylab="Pathway expression level", bty="l", cex.main = 1, xaxt = "n", shade= F, se=3, scheme=1)
      xmin <- par("usr")[1]
      xmax <- par("usr")[2]
      points(x=plot_df$time, y=plot_df$value, pch=20, col="navy", cex=0.9, xpd=NA)
      text(x=plot_df$time, y=plot_df$value, labels=plot_df[,1], pos = 4, cex = 0.5, col="navy", xpd=NA)
      legend("topright", legend = "Cluster", pch = 20, col = "navy", bty="n", text.col="navy", cex=0.9)
      legend("topright", legend=paste0("\nAdjusted p-value = ", round(varying_pathways[[i]], 5)), bty="n", cex=0.9)
      axis(side=1, at=c(xmin, xmax), labels = c("Early", "Late"), tick=T)


    dev.off()
      
  }

}                         
                         
IdentifyVaryingPWs <- function(object, pval_threshold=0.05){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@n.pcs)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  if (is.null(object@cluster.pathways)){
    stop("CalculatePWProfiles has not been run. See ?Tempora::CalculatePWProfiles for details")
  }
  gsva_bycluster <- object@cluster.pathways

  significant_pathways <- c()
  for (i in 1:object@n.pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways <- c(names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)), significant_pathways)
  }

  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways <- gsub("\\s*\\([^\\)]+\\)","",pca_pathways)
  pca_pathways_cleaned <- gsub("[[:punct:]]", "", pca_pathways)
  themes <- pca_pathways_cleaned

  cat("Fitting GAM models...")

  p_vals <- gams <- list()
  for (i in 1:length(themes)){
      if (i %% 1000 == 0)
      {
        print(paste(i, "of", length(themes)))
      }
    
    if(length(grep(themes[i], rownames(gsva_bycluster))) == 0) {
      p_vals[[i]] <- 1
      gams[[i]] <- NA
      next
      }
    if (length(grep(themes[i], rownames(gsva_bycluster))) > 1){
      plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ], na.rm=T))
    } else if (length(grep(themes[i], rownames(gsva_bycluster))) == 1){
      plot_df <- data.frame(cluster=names(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]) }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams[[i]] <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams[[i]])
    p_vals[[i]] <- temp_anova$s.pv
  }

  names(p_vals) <- names(gams) <- themes

  pval_threshold = pval_threshold
  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]

  if (length(varying_pathways)==0){
    cat("No temporally varying pathways detected. Please try running IdentifyVaryingPWs with a more relaxed p-value cutoff.")
    #eventhough the function was not successful return the object because in the vignette
    # this function call sets the original object to what is returned and if it is null
    # you loose all the processing you have done until now. 
    return(object)
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
}


makeTrajectoryPlot = function(covid_tempora, colors=NULL)
{
  edge_graph <- igraph::graph_from_data_frame(d=covid_tempora@trajectory, vertices = covid_tempora@cluster.metadata, directed = T) 
l <- igraph::layout_with_sugiyama(edge_graph, hgap=50, vgap=100, layers = covid_tempora@cluster.metadata$Cluster_time_score, maxiter = 2000)
colours <- brewer.pal(length(levels(covid_tempora@meta.data$Timepoints)), "YlOrRd")

g = edge_graph
# precompute layout
V(g)$x <- l$layout[,1]
V(g)$y <- l$layout[,2]
    
    print(l)
    
V(g)$label=unlist(lapply(str_split(covid_tempora@cluster.metadata$label, "-"), function(x){x[2]}))
aesDF=data.frame("x"=V(g)$x, "y"=V(g)$y)
for (tpx in unique(covid_tempora@meta.data$Timepoints))
{
  tpData = unlist(lapply(1:nrow(covid_tempora@cluster.metadata), function(x) as.numeric(covid_tempora@cluster.metadata[x,2:((length(levels(covid_tempora@meta.data$Timepoints)))+1)][tpx])))
  aesDF[[as.character(tpx)]] = tpData
}

xMin = round(min(l$layout[,1]))
xMax = round(max(l$layout[,1]))

xMin = xMin - 0.05*(xMax-xMin)
xMax = xMax + 0.05*(xMax-xMin)

yMin = round(min(l$layout[,2]))
yMax = round(max(l$layout[,2]))

yMin = yMin - max(c(0.05*(yMax-yMin), 100))
yMax = yMax + max(c(0.05*(yMax-yMin), 100))

print(paste(xMin, xMax, yMin, yMax))

if (is.null(colors))
{
getPalette = colorRampPalette(brewer.pal(10, "Set3"))
nTimepoints = length(unique(covid_tempora@meta.data$Timepoints))
colors=getPalette(nTimepoints)
}                         


aesCols = as.character(levels(unique(covid_tempora@meta.data$Timepoints)))

#geom_scatterpie(cols = aesCols, data = aesDF, colour = "white", pie_scale = 0.5) +
                                                  
p=ggraph(g, "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link(arrow = arrow(length = unit(5, 'mm')), end_cap = circle(7, 'mm')) + xlim(xMin, xMax)+ geom_scatterpie(cols = aesCols, data = aesDF, colour = "white", pie_scale = 0.5) +
  geom_text_repel(aes(x=x, y=y, label=label), point.size=20, size=5, max.overlaps=10)+
  scale_fill_manual(name="Timepoint", values = colors) +
  coord_fixed() + 
  scale_y_continuous(expand = c(0, 0),limits = c(yMin, yMax), breaks = c(yMin,yMax), minor_breaks = NULL, label=c("Late", "Early"))+ ylab("Inferred Time")+
  theme_minimal() +  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face="bold", size=18),
    axis.title.y = element_text(face="bold", size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color='black'),
    legend.text=element_text(size=24),
    legend.title = element_text(size=24)
    )

return(p)
}

p=makeTrajectoryPlot(tempora.all.covid, colors=c("#CCCFCF", "#F6DFEE","#E8AFD1","#EE7DAF"))
save_plot(p, paste("stateCellnamesTP_all_covid", "trajectory", sep="_"),  fig.width = 30, fig.height = 10)



makeTrajectoryPlot = function(covid_tempora, colors=NULL)
{
  edge_graph <- igraph::graph_from_data_frame(d=covid_tempora@trajectory, vertices = covid_tempora@cluster.metadata, directed = T) 
l <- igraph::layout_with_sugiyama(edge_graph, hgap=50, vgap=100, layers = covid_tempora@cluster.metadata$Cluster_time_score, maxiter = 2000)
colours <- brewer.pal(length(levels(covid_tempora@meta.data$Timepoints)), "YlOrRd")

g = edge_graph
# precompute layout
V(g)$x <- l$layout[,1]
V(g)$y <- l$layout[,2]
    
    
                     
l$layout[,1] = 500*(l$layout[,1]/max(l$layout[,1]))

xMinO = round(min(l$layout[,1]))
xMaxO = round(max(l$layout[,1]))

xMin = xMinO - max(c(0.05*(xMaxO-xMinO), 100))
xMax = xMaxO + max(c(0.05*(xMaxO-xMinO), 100))

yMinO = round(min(l$layout[,2]))
yMaxO = round(max(l$layout[,2]))

yMin = yMinO - max(c(0.05*(yMaxO-yMinO), 100))
yMax = yMaxO + max(c(0.05*(yMaxO-yMinO), 100))
                         
V(g)$x <- l$layout[,1]
V(g)$y <- l$layout[,2]
    
V(g)$label=unlist(lapply(str_split(covid_tempora@cluster.metadata$label, "-"), function(x){x[2]}))
aesDF=data.frame("x"=V(g)$x, "y"=V(g)$y)
for (tpx in unique(covid_tempora@meta.data$Timepoints))
{
  tpData = unlist(lapply(1:nrow(covid_tempora@cluster.metadata), function(x) as.numeric(covid_tempora@cluster.metadata[x,2:((length(levels(covid_tempora@meta.data$Timepoints)))+1)][tpx])))
  aesDF[[as.character(tpx)]] = tpData
}
aesDF$radius = 25

print(l)

print(paste(xMin, xMax, yMin, yMax))

if (is.null(colors))
{
getPalette = colorRampPalette(brewer.pal(10, "Set3"))
nTimepoints = length(unique(covid_tempora@meta.data$Timepoints))
colors=getPalette(nTimepoints)
}                         


aesCols = as.character(levels(unique(covid_tempora@meta.data$Timepoints)))

#geom_scatterpie(cols = aesCols, data = aesDF, colour = "white", pie_scale = 0.5) +
                                                  
p=ggraph(g, "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link(arrow = arrow(length = unit(5, 'mm')), end_cap = circle(7, 'mm')) + xlim(xMin, xMax)+geom_scatterpie(aes(x=x, y=y, r=radius), cols = aesCols, data = aesDF, colour = "white") +
  geom_text_repel(aes(x=x, y=y, label=label), point.size=20, size=5, max.overlaps=10)+
  scale_fill_manual(name="Timepoint", values = colors) +
  coord_fixed() + 
  scale_y_continuous(expand = c(0, 0),limits = c(yMin, yMax), breaks = c(yMin,yMax), minor_breaks = NULL, label=c("Late", "Early"))+ ylab("Inferred Time")+
  theme_minimal() +  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face="bold", size=18),
    axis.title.y = element_text(face="bold", size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color='black'),
    legend.text=element_text(size=24),
    legend.title = element_text(size=24)
    )

return(p)
}


analyseTempora = function( scobj, clusterField, clusterFieldLevels, timepointField, timepoints, outprefix, pw.threshold=0.05, pw.max=50, pw.genes.max=500,pw.genes.min=10,do.pw=TRUE,colors=NULL )
{
  covid_tempora <- ImportSeuratObject(scobj,
                                     clusters = clusterField,
                                     timepoints = timepointField, 
                                     assayType = "RNA",
                                     cluster_labels=clusterFieldLevels,
                                     timepoint_order = timepoints)

gmt_file = "Human_COVID.gmt"


covid_tempora <- CalculatePWProfiles(covid_tempora, 
                gmt_path = gmt_file,
                method="gsva", min.sz = pw.genes.min, max.sz = pw.genes.max, parallel.sz = 1)
n_calc_pcs = ncol(covid_tempora@cluster.pathways.dr$x)
print(paste("Calculated number of PCs:", n_calc_pcs))


print("Building Trajectory")
covid_tempora = Tempora::BuildTrajectory(covid_tempora, n_pcs = n_calc_pcs, difference_threshold = 0.01)



print("Drawing Trajectory")

p = makeTrajectoryPlot(covid_tempora,colors)
p=save_plot(p, paste(outprefix, "trajectory", sep="_"),  fig.width = 15, fig.height = 5)

if (do.pw)
{
print("IdentifyVaryingPWs")
#Fit GAMs on pathway enrichment profile
covid_tempora <- IdentifyVaryingPWs(covid_tempora, pval_threshold = pw.threshold)

print("PlotVaryingPWs")
PlotVaryingPWs(covid_tempora, clusterField, outprefix, pw.max)    
}


return(covid_tempora)
}


makeStateCellnamesTP = function( scobj, thresh=3 )
{
  stateCellnamesTPnp = paste(unlist(
    lapply(str_split(scobj$unified_cellnames, ";"), function(x){return(str_replace(str_replace(x[1], "Gamma delta T cells", "gd T cells"),
     "Plasmacytoid dendritic cell", "pDC"))})
  ), substr(scobj$disease_state, 0,4), scobj$ctrlTPs)

  print(table(stateCellnamesTPnp))

  scobj$stateCellnamesTP = stateCellnamesTPnp

  cellThreshold = thresh
  keepClusters = names(table(stateCellnamesTPnp)[table(stateCellnamesTPnp) > cellThreshold])
  removeClusters = table(stateCellnamesTPnp)[table(stateCellnamesTPnp) <= cellThreshold]
  print("removing clusters")
  print(removeClusters)

  scobj = subset(scobj,stateCellnamesTP %in% keepClusters)

  stateCellnamesTPnp = paste(unlist(
    lapply(str_split(scobj$unified_cellnames, ";"), function(x){return(str_replace(str_replace(x[1], "Gamma delta T cells", "gd T cells"),
     "Plasmacytoid dendritic cell", "pDC"))})
  ), substr(scobj$disease_state, 0,4), scobj$ctrlTPs)

  stateCellnamesTPFactornp = factor(stateCellnamesTPnp)
  stateCellnamesTPNumericnp = as.numeric(stateCellnamesTPFactornp)
  scobj$stateCellnamesTP = stateCellnamesTPNumericnp

  return(list("seurat"=scobj, "factorObj"=stateCellnamesTPFactornp))
}


#Load MouseCortex sample data
obj.integrated = readRDS("seurat_object_for_tempora.Rds")

ctrlTPs = obj.integrated$mpoint
ctrlTPs[obj.integrated$disease_state == "CONTROL"] = "0"
obj.integrated$ctrlTPs = ctrlTPs
unique(ctrlTPs)

# create new (numeric) clusters for tempora
# disease_state => unique(obj.integrated$disease_state)
# unified_cellnames => as.character(unique(obj.integrated$unified_cellnames))

#niceCellnames = unlist(
  #lapply(str_split(obj.integrated$unified_cellnames, ";"), function(x){return(x[1])})
#)



#stateCellnames = paste(niceCellnames, substr(obj.integrated$disease_state, 0,4))
#stateCellnamesFactor = factor(stateCellnames)
#stateCellnamesNumeric = as.numeric(stateCellnamesFactor)

#obj.integrated$stateCellnames = stateCellnamesNumeric


#stateCellnamesTP = paste(niceCellnames, substr(obj.integrated$disease_state, 0,4), obj.integrated$ctrlTPs)
#stateCellnamesTPFactor = factor(stateCellnamesTP)
#stateCellnamesTPNumeric = as.numeric(stateCellnamesTPFactor)
#obj.integrated$stateCellnamesTP = stateCellnamesTPNumeric


procObj.symp_only = makeStateCellnamesTP(subset(obj.integrated, disease_state=="SYMPTOMATIC"))
tempora.symp_only = analyseTempora(procObj.symp_only$seurat, "stateCellnamesTP", levels(procObj.symp_only$factorObj), "ctrlTPs", c("1", "2", "3"), "stateCellnamesTP_symp_only", pw.threshold=0.5)

procObj.asymp_only = makeStateCellnamesTP(subset(obj.integrated, disease_state=="ASYMPTOMATIC"))
tempora.asymp_only = analyseTempora(procObj.asymp_only$seurat, "stateCellnamesTP", levels(procObj.asymp_only$factorObj), "ctrlTPs", c("1", "2", "3"), "stateCellnamesTP_asymp_only", pw.threshold=0.5)


procObj.asymp_symp_only = makeStateCellnamesTP(subset(obj.integrated, disease_state!="CONTROL"))
tempora.asymp_symp_only = analyseTempora(procObj.asymp_symp_only$seurat, "stateCellnamesTP", levels(procObj.asymp_symp_only$factorObj), "ctrlTPs", c("1", "2", "3"), "stateCellnamesTP_asymp_symp_only", pw.threshold=0.5)



procObj.all = makeStateCellnamesTP(obj.integrated)
tempora.all = analyseTempora(procObj.all$seurat, "stateCellnamesTP", levels(procObj.all$factorObj), "ctrlTPs", c("0", "1", "2", "3"), "stateCellnamesTP_all", pw.threshold=0.5)

procObj.symp = makeStateCellnamesTP(subset(obj.integrated, disease_state=="SYMPTOMATIC" | disease_state=="CONTROL"))
tempora.symp = analyseTempora(procObj.symp$seurat, "stateCellnamesTP", levels(procObj.symp$factorObj), "ctrlTPs", c("0", "1", "2", "3"), "stateCellnamesTP_symp", pw.threshold=0.5)

procObj.asymp = makeStateCellnamesTP(subset(obj.integrated, disease_state=="ASYMPTOMATIC" | disease_state=="CONTROL"))
tempora.asymp = analyseTempora(procObj.asymp$seurat, "stateCellnamesTP", levels(procObj.asymp$factorObj), "ctrlTPs", c("0", "1", "2", "3"), "stateCellnamesTP_asymp", pw.threshold=0.5)


procObj.all.covid = makeStateCellnamesTP(obj.integrated)
tempora.all.covid = analyseTempora(procObj.all.covid$seurat, "stateCellnamesTP",
levels(procObj.all.covid$factorObj), "ctrlTPs", c("0", "1", "2", "3"), "stateCellnamesTP_all_covid", pw.threshold=0.8)

PlotVaryingPWs(tempora.all.covid, "stateCellnamesTP","stateCellnamesTP_all_covid")
p=makeTrajectoryPlot(tempora.all.covid)
save_plot(p, paste("stateCellnamesTP_all_covid", "trajectory", sep="_"),  fig.width = 24, fig.height = 16)



makeStateCellnames = function( scobj, thresh=3 )
{
  stateCellnamesTPnp = paste(unlist(
    lapply(str_split(scobj$unified_cellnames, ";"), function(x){return(x[1])})
  ), substr(scobj$disease_state, 0,4))

  print(table(stateCellnamesTPnp))

  scobj$stateCellnamesTP = stateCellnamesTPnp

  cellThreshold = thresh
  keepClusters = names(table(stateCellnamesTPnp)[table(stateCellnamesTPnp) > cellThreshold])
  removeClusters = table(stateCellnamesTPnp)[table(stateCellnamesTPnp) <= cellThreshold]
  print("removing clusters")
  print(removeClusters)

  scobj = subset(scobj,stateCellnamesTP %in% keepClusters)

  stateCellnamesTPnp = paste(unlist(
    lapply(str_split(scobj$unified_cellnames, ";"), function(x){return(x[1])})
  ), substr(scobj$disease_state, 0,4))
  stateCellnamesTPFactornp = factor(stateCellnamesTPnp)
  stateCellnamesTPNumericnp = as.numeric(stateCellnamesTPFactornp)
  scobj$stateCellnames = stateCellnamesTPNumericnp

  return(list("seurat"=scobj, "factorObj"=stateCellnamesTPFactornp))
}

procObj.allCT = makeStateCellnames(obj.integrated)
tempora.allCT = analyseTempora(procObj.allCT$seurat, "stateCellnames", levels(procObj.allCT$factorObj), "ctrlTPs", c("0", "1", "2", "3"), "stateCellnames_all", pw.threshold=0.5)


p=makeTrajectoryPlot(tempora.allCT, colors=c("#CCCFCF", "#F6DFEE","#E8AFD1","#EE7DAF"))
save_plot(p, paste("stateCellnames_all", "trajectory", sep="_"),  fig.width = 30, fig.height = 10)


save.image("covid_tempora.RData")
load("covid_tempora.RData")



selCelltypes = c("Monocytes", "B cells", "CD4+ T cells", "CD8+ T cells", "NK cells")


makecelltypenamesTP = function( scobj, thresh=3 )
{
  stateCellnamesTPnp = paste(scobj$celltypenames, substr(scobj$disease_state, 0,4), scobj$ctrlTPs)

  print(table(stateCellnamesTPnp))

  scobj$stateCellnamesTP = stateCellnamesTPnp

  cellThreshold = thresh
  keepClusters = names(table(stateCellnamesTPnp)[table(stateCellnamesTPnp) > cellThreshold])
  removeClusters = table(stateCellnamesTPnp)[table(stateCellnamesTPnp) <= cellThreshold]
  print("removing clusters")
  print(removeClusters)

  scobj = subset(scobj,stateCellnamesTP %in% keepClusters)

  stateCellnamesTPnp = paste(scobj$celltypenames, substr(scobj$disease_state, 0,4), scobj$ctrlTPs)

  stateCellnamesTPFactornp = factor(stateCellnamesTPnp)
  stateCellnamesTPNumericnp = as.numeric(stateCellnamesTPFactornp)
  scobj$stateCellnamesTP = stateCellnamesTPNumericnp

  return(list("seurat"=scobj, "factorObj"=stateCellnamesTPFactornp))
}



for (selCelltype in selCelltypes)
{


procObj.symp.selCelltype = makecelltypenamesTP(subset(obj.integrated, celltypenames %in% selCelltype))
tempora.symp.selCelltype = analyseTempora(procObj.symp.selCelltype$seurat, "stateCellnamesTP",
levels(procObj.symp.selCelltype$factorObj), "ctrlTPs", c("0", "1", "2", "3"), paste("tempora_sc/stateCellnamesTP_all_", selCelltype, sep=""), pw.threshold=0.8, pw.max=50, pw.genes.max=250, do.pw=T,colors=c("#CCCFCF", "#F6DFEE","#E8AFD1","#EE7DAF"))

procObj.symp.selCelltype = makecelltypenamesTP(subset(obj.integrated, celltypenames %in% selCelltype & (disease_state=="SYMPTOMATIC" | disease_state=="CONTROL")))
tempora.symp.selCelltype = analyseTempora(procObj.symp.selCelltype$seurat, "stateCellnamesTP",
levels(procObj.symp.selCelltype$factorObj), "ctrlTPs", c("0", "1", "2", "3"), paste("tempora_sc/stateCellnamesTP_sympt_", selCelltype, sep=""), pw.threshold=0.8, pw.max=50, pw.genes.max=250, do.pw=T,colors=c("#CCCFCF", "#F6DFEE","#E8AFD1","#EE7DAF"))

procObj.symp.selCelltype = makecelltypenamesTP(subset(obj.integrated, celltypenames %in% selCelltype & (disease_state=="ASYMPTOMATIC" | disease_state=="CONTROL")))
tempora.symp.selCelltype = analyseTempora(procObj.symp.selCelltype$seurat, "stateCellnamesTP",
levels(procObj.symp.selCelltype$factorObj), "ctrlTPs", c("0", "1", "2", "3"), paste("tempora_sc/stateCellnamesTP_asympt_", selCelltype, sep=""), pw.threshold=0.8, pw.max=50, pw.genes.max=250, do.pw=T,colors=c("#CCCFCF", "#F6DFEE","#E8AFD1","#EE7DAF"))


}


library(Seurat)
library(writexl)
library(data.table)
library(stringr)

obj.integrated = readRDS("seurat_object_for_tempora.Rds")


makesummary = function(a, suffix)
{
  out = {}
  out["num"] = length(a)
  
  if (length(a) == 0)
  {
    f = c(0,0,0,0,0)
    meanA = 0
    sumA = 0
  } else {
    f = fivenum(a)
    meanA = mean(a)
    sumA = sum(a)
  }
  out["min"] = f[1]
  out["lower_hinge"] = f[2]
  out["median"] = f[3]
  out["upper_hinge"] = f[4]
  out["max"] = f[5]
  out["mean"] = meanA
  out["sum"] = sumA
  
  names(out) = paste(names(out), suffix, sep=".")
  
  return(out)
}


getExprDataPop = function(markerObj, markerCells, sampleSuffix, slot="data")
{
  expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = slot)
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

get_population_expression_data = function(scobj, group, outname)
{
  
  exprData = list()
  for (cellPop in unique(as.character(unlist(scobj[[group]]))))
  {
    print(cellPop)
    varPop = str_to_lower( str_replace_all(
                            str_replace_all(#
                              str_replace_all( cellPop, "\\(|\\)| |,", "_"),
                              "__", "_"),
                            "_$", "")
                           )
    
    
    cellPopCells = rownames(scobj[[group]][scobj[[group]] == cellPop])
  
    exprData[[varPop]] = getExprDataPop(markerObj=scobj, markerCells=cellPopCells, sampleSuffix=varPop, slot="counts")
  }
  
  
  meanExprData = list()
  
  for (name in names(exprData))
  {
    
    exprDf = as.data.frame(exprData[[name]])
    subdf = exprDf[ ,c("gene", paste("sum", name, sep=".")) ]
  
    meanExprData[[name]] = subdf
  }
  
  cellnames_manualExprDF = Reduce(function(x,y) merge(x = x, y = y, by = "gene"), meanExprData)
  
  write.table(cellnames_manualExprDF, file = paste(outname, ".tsv", sep=""), quote=FALSE, sep = "\t", row.names = F)
  write_xlsx( cellnames_manualExprDF, path = paste(outname, ".xlsx", sep="") )
  
  return(cellnames_manualExprDF)
  
  
}


manualExprDF = get_population_expression_data(obj.integrated, "celltypenamestp", "population_expr/sum_exprdf_celltypenamestp")

manualExprDF = get_population_expression_data(subset(obj.integrated, cells=names(obj.integrated$disease_state[obj.integrated$disease_state=="ASYMPTOMATIC"])), "celltypenamestp", "population_expr/sum_exprdf_celltypenamestp_asympt")
manualExprDF = get_population_expression_data(subset(obj.integrated, cells=names(obj.integrated$disease_state[obj.integrated$disease_state=="SYMPTOMATIC"])), "celltypenamestp", "population_expr/sum_exprdf_celltypenamestp_sympt")
manualExprDF = get_population_expression_data(subset(obj.integrated, cells=names(obj.integrated$disease_state[obj.integrated$disease_state=="CONTROL"])), "celltypenamestp", "population_expr/sum_exprdf_celltypenamestp_ctrl")



save.image("covid_tempora_end.RData")
