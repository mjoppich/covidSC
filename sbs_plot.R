
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

make_descr_label = function(plot, descr)
{
  descrLabel <- ggdraw() + draw_label(descr, fontface='bold', angle = 0)
  
  pe = cowplot::plot_grid(descrLabel, plot, ncol=1, nrow=2, labels=NULL,rel_heights = c(0.1, 1),
                          align = "h", axis = "l")
  
  return(pe)
}

makeSideBySideDotPlot = function(scobj, plotElems, featureGenes = c("Csf1r", "Chil3"), group.by="cellnames_manual",   col.min = -3, col.max = 3, cols = c("grey", "blue"), scaled=T, features.rotate=T, rel.width.legend=0.05, plot.title="")
{
  
  featureGenes.orig = featureGenes
  featureGenes = intersect(featureGenes, rownames(scobj))
  
  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  allFoundGroups = c()
  
  elemCount = 0
  for (plotName in names(plotElems))
  {
    print(plotName)
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    if (scaled)
    {
      plotElem_orig = DotPlot(subset(scobj, cells=plotCells), features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols=cols, scale=scaled)
      
    } else {
      plotElem_orig = DotPlot(subset(scobj, cells=plotCells), features=featureGenes, group.by = group.by, cols=cols, scale=scaled)
      
      # https://github.com/satijalab/seurat/issues/2798
      plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
      
      plotElem_orig$data[plotElem_orig$data$avg.exp>col.max, 'avg.exp'] = col.max
      plotElem_orig$data[plotElem_orig$data$avg.exp<col.min, 'avg.exp'] = col.min
      
      plotElem_orig <- plotElem_orig + 
        geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
        guides(color = guide_colorbar(title = 'Average Expression'))
    }
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    allFoundGroups = unique(c(allFoundGroups, as.character(plotElem_orig$data$id)))
    
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    if (features.rotate)
    {
      plotElem_orig= plotElem_orig + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
    }
    
    plot_orig_list[[plotName]] = plotElem_orig
  }
  
  
  if (scaled)
  {
    
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
      
      allUnscaledValues = data.frame()
      for (plotName in names(plot_orig_list))
      {
        pData = plot_orig_list[[plotName]]$data
        newDF = data.frame(pData[ pData$features.plot==featureName, ]$avg.exp)
        colnames(newDF) = c(plotName)
        
        allUnscaledValues = cbind.fill(allUnscaledValues, newDF)
      }
      
      sampleNames = colnames(allUnscaledValues)
      allUnscaledValues = as.data.frame(allUnscaledValues)
      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, sampleNames)
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      allScaledValues
      
      for (plotName in names(plot_orig_list))
      {
        
        plotElem_orig = plot_orig_list[[plotName]]
        pData = plotElem_orig$data
        
        
        # https://github.com/satijalab/seurat/issues/2798
        plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
        
        plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp.scaled2"] = allScaledValues[!is.na(allScaledValues[, plotName]),plotName]
        plot_orig_list[[plotName]] = plotElem_orig
      }
    }
    
    
    
    for (plotName in names(plot_orig_list))
    {
      for (featureName in allFoundFeatures)
      {
        plotElem_orig = plot_orig_list[[plotName]]
        pData = plotElem_orig$data
        
        for (groupName in allFoundGroups)
        {
          
          subdf = pData[pData$id == groupName & pData$features.plot == featureName,]
          
          if (nrow(subdf) == 0)
          {
            print(paste(plotName, featureName, groupName))
            addDF = data.frame(avg.exp=0, pct.exp=0, features.plot=featureName, id=groupName, avg.exp.scaled=col.min, avg.exp.scaled2=col.min, stringsAsFactors = FALSE)
            rownames(addDF) = featureName
            
            plotElem_orig$data$id <- as.character(plotElem_orig$data$id)
            
            plotElem_orig$data = rbind(plotElem_orig$data, addDF)
            
            plotElem_orig$data$id <- as.factor(plotElem_orig$data$id)
            
            plot_orig_list[[plotName]] = plotElem_orig
          }
          
        }
      }
    }
    
    for (plotName in names(plot_orig_list))
    {
      plotElem_orig = plot_orig_list[[plotName]]
      pData = plotElem_orig$data
      
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min
      
      plotElem_orig <- plotElem_orig + 
        geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled2')) +
        guides(color = guide_colorbar(title = 'Average Scaled Expression'))
      
      plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
      plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
      
      plot_orig_list[[plotName]] = plotElem_orig
    }
  }
  
  
  
  elemCount = 0
  plot_list=list()
  for (plotName in names(plot_orig_list))
  {
    
    plotElem_orig = plot_orig_list[[plotName]]
    
    elemCount = elemCount + 1
    
    plotElem = plotElem_orig +  theme(
      axis.text.y = element_blank(),
      panel.grid.major.x = element_line( size=.1, color="black" ),
      axis.line.y = element_line(size = 0),
      axis.ticks.length.y = unit(0, "points")
    )
    
    pe = make_descr_label(plotElem, plotElems[[plotName]]$label)
    
    plot_list[[plotName]] = pe
    
  }
  
  legendDescr = 'Average Expression'
  if (scaled)
  {
    legendDescr = 'Average Scaled Expression'  
  }
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    plotElem_orig + 
      guides(color = guide_colorbar(title = legendDescr, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )
  
  title <- ggdraw() + draw_label(plot.title, fontface='bold')
  
  plotElem_new = plot_orig_list[[names(plot_orig_list)[[1]]]]
  plotElem_new$layers[[1]] = NULL
  pa = plotElem_new + theme(axis.text.x=element_text(colour = "white")) #+ xlim("", "B")
  pal = make_descr_label(pa, "")
  
  
  plist = list()
  plist[["Legend"]] = pal
  
  for (le in names(plot_list))
  {
    plist[[le]] = plot_list[[le]]
  }
  
  ap=cowplot::plot_grid(
    plotlist = plist,
    labels = NULL,
    nrow=1, rel_widths = c(rel.width.legend, rep(0.1, length(plot_list))),
    align = "v", axis="bt"
  )
  
  fp = cowplot::plot_grid(title, ap, legend_b, ncol = 1, rel_heights = c(0.05, 1, 0.1) )
  return(fp)
}