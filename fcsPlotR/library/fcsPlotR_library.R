library(data.table)
library(ggplot2)
library(reshape2)
library(tools)
library(flowCore)
library(plyr)
library(doMC)
library(boot)
registerDoMC()
getDoParWorkers()

#### Helper functions for loading FCS files ####


#' Wrapper for flow Core read.FCS
#'
#' @param filePath path to the FCS file
#' @return A flowframe object
#' @examples
#' read.FCS('/test/testfcs.fcs')
bb.loadFCS <-function(filePath, ...){
  # a wrapper to the standard read.FCS, optimized for cyTOF data 
  datFCS <- read.FCS(filePath,min.limit=NULL,transformation = 'linearize', ...)  
  return(datFCS)
}

#' Converts FlowFrame to data table
#' 
#' @param datFCS A flowframe e.g. from read.FCS or bb.loadFCS
#' @return A data.table with row=cell, column=Channel
#' @examples
#' read.FCS('/test/testfcs.fcs')
bb.flowFrame2dt <-function(datFCS){
  # converts a flow frame from read.FCS to a data table
  dt = datFCS@exprs
  colnames(dt) = make.names(make.unique(colnames(dt)))
  fil = !is.na(datFCS@parameters@data$desc)
  colnames(dt[,fil]) <- datFCS@parameters@data$desc[fil]
  dt <- data.table(dt)
  uninames = make.names(make.unique(datFCS@parameters@data$name))
  uni_newnames = make.names(make.unique(datFCS@parameters@data$desc))
  setnames(dt,uninames[fil],uni_newnames[fil])
  return(dt)
}

#' Loads multiple FCS files at once
#'
#' @param fileList List of files to be loaded
#' @param fileDir the base directory for the files or a directory with FCS files
#' @param condDict a named vector which maps Filenames to condition
#' 
#' @return a dataframe with all combined FCS files
bb.loadConvertMultiFCS <- function(fileList=NaN,fileDir=NaN,condDict=NaN,subSample = NA){
  if (is.na(fileList)){
    fileList = list.files(fileDir, pattern='*.fcs')
  }
  dat = data.table()
  for(file in fileList){
    if(file_ext(file) =='fcs'){
      if (!is.na(condDict)){
        cond = condDict[[file]]
      } else {
        cond = file
      }
      
      #load and convert to dt
      tmpDT <- bb.flowFrame2dt(bb.loadFCS(file.path(fileDir,file)))
      
      if (!is.na(subSample) && is.numeric(subSample)){
        if (subSample < nrow(tmpDT)){
          tmpDT = tmpDT[sample(1:nrow(tmpDT),subSample,replace = F), ]
        } else {stop('Less cells than subSample size!')}
      }      
      
      #combine dt
      dat = rbindlist(list(dat,tmpDT[,'condition':=cond]), fill=T)
    }
  }
  return(dat)
}

#' Generates a metadata table from a
#'
#'
#bb.flowFrame2metaDt <-function(datFCS){
#  # gets the metadata table
#  metaDt <- data.table(datFCS@parameters@data)
#  return(metaDt)
#}


#' Extracts information from string fields
#' 
#' @param name a string
#' @param seq the seperator that seperates fields in the string
#' @param strPos numeric or vector: indicates which positions of the string should be extracted
#' @param censorStr string that gets ignored
#' @return The extracted stringfield
#'
bb.getInfoFromString<- function(name,sep='_',strPos=2,censorStr='.fcs'){
  tp <- gsub(censorStr,'',strsplit(name,sep)[[1]][strPos])
  tp <- paste(tp,collapse = sep)
  return(tp)
}


#' Parses a list of strings and gets the corresponding information
#' 
#' See \code(bb.getInfoFromString) for more information
#'
bb.getInfoFromFileList <- function(fileList,sep='_',strPos=2,censorStr='.fcs'){
  condDict = sapply(fileList,function(file){bb.getInfoFromString(file,sep,strPos,censorStr)})
  names(condDict) <- fileList
  return(condDict)
}



#' Calculates the bhSNE from a melted data table
#'
#' @param input_dat melted data_table
#' @param channels list of channels to be Used
#' @param channel_var column with the channel names
#' @param value_var column to be used as values
#' @param id_var column name for the unique cell ID
#' @param group_var variable that groups the table
#' @param scale should the input data be rescaled? 
#' @param subsample_groups T: resamples all groups to the group with the least members, Integer: resamples all groups to the integer
#' @param subsample_mode 
#'  'equal': if integer provided to resample_groups all groups will contain min(resample_groups, min_groupsize) cells
#'  'unequal': if integer provided to resample_groups, all groups will contain maximally resample_groups cells (or less)
#' @param verbose: should bhSne be verbose?
#' 
#' @return 
bb.calcTSNE <- function(input_dat, channels, value_var='counts', channel_var='channel',
                        id_var='id', group_var ='condition', scale=F,
                        subsample_groups=F, subsample_mode='equal',
                        verbose=T,
                        dims=2,
                        ...){

  # do subsampling
  if (is.numeric(subsample_groups)){
    if (subsample_mode == 'equal'){
      min_n = min(input_dat[get(channel_var) == channels[1], .(n= .N), by=get(group_var)]$n)
      subsample_groups = min(subsample_groups, min_n)
    }
    ids = input_dat[get(channel_var) == channels[1], .(fil = get(id_var)[sample.int(.N, min(subsample_groups, .N))]),by=get(group_var)]$fil
    
  } else if (subsample_groups){
    min_n = min(input_dat[get(channel_var) == channels[1],.(n= .N), by=get(group_var)]$n)
    ids = input_dat[get(channel_var) == channels[1], .(fil = get(id_var)[sample.int(.N, min_n)]),by=get(group_var)]$fil
    
  } else {
    ids = input_dat[, get(id_var)]
  }


  
  dt = dcast.data.table(input_dat[(get(channel_var) %in% good_channels) & get(id_var) %in% ids], formula = as.formula(paste(id_var, '~', channel_var)), value.var = value_var)
  
  if (scale){
    tsnedat = scale(dt[, channels, with=F])
    
  } else {
    tsnedat = dt[, channels, with=F]
  }
  
  tsne_out <- Rtsne(tsnedat, verbose=verbose,dims=dims,...)
  tsne_out$Y = data.table(tsne_out$Y)
  setnames(tsne_out$Y, names(tsne_out$Y), paste("bh",names(tsne_out$Y),sep='_'))
  tsne_out$Y$id = dt[, get(id_var)]
  setnames(tsne_out$Y, 'id', id_var)
  setkeyv(tsne_out$Y, id_var)
  tsne_out$channels = channels
  tsne_out$scale= F
  if (group_var %in% names(input_dat)){
    tsne_out$groups= input_dat[, unique(get(group_var))]
  } else {
    tsne_out$groups = NaN
  }
  
  tsne_out$subsample_groups = subsample_groups
  tsne_out$subsample_mode = subsample_mode
  return(tsne_out)

}

### Network plotting functions ####

#' Loads a Cytoscape JSON as an igraph graph
#' 
#' @param  jdat a loaded Cytoscape JSON
#' @result an igraph-graph
#' 
#' @example 
#' jdat = fromJSON(file = 'xy.json')
#' g = load_json_graph(jdat)
#' list.edge.attributes(g)
#' list.vertex.attributes(g)
#' 
#' plot.igraph(g,
#'            layout=-cbind(V(g)$posx, V(g)$posy),
#'            rescale=T,
#'            edge.arrow.size=0.1,
#'            vertex.size =7,
#'            vertex.label.cex = 0.5)
bb.load_json_graph = function(jdat){
  "
  Loads a graph with metadata from a cytoscape.js JSON file
  "
  edgelist = sapply(jdat$elements$edges, function(e){
    # Loads the edge matrix f  
    return(c(e$data$source, e$data$target))
  })
  # Create graph
  g = graph_from_edgelist(t(edgelist))
  
  # Loop over the nodes to get metadata
  for (v in jdat$elements$nodes){
    
    idx = which(V(g)$name == v$data$id)
    
    # Loop over the metadata entries
    for (at in names(v$data)){
      
      g = set.vertex.attribute(g, at, idx, v$data[[at]])
    }
    
    # set position metadata
    g = set.vertex.attribute(g,'posx',idx, as.numeric(v$position$x))
    g = set.vertex.attribute(g,'posy',idx, as.numeric(v$position$y))
    
  }
  
  # Loop over edges metadata
  for (idx in 1:length(jdat$elements$edges)){
    e = jdat$elements$edges[[idx]]
    for (at in names(e$data)){
      g = set.edge.attribute(g,at,idx, e$data[[at]])
    }
  }
  return (g)
}





### various small helper functions ####

#' Maps a vector x to a colormap and returns the RGB values
#'
#' @param x values to be mapped
#' @param xmax maximum value to be mapped
#' @param symmetric if T the map will be symmetric from -xmax to xmax centered arround 0
#' @param cmap a colormap, e.g. from colorRamp
#' @param na_col what color should NA have (name or RGB)
#' 
#' @return returns the values mapped
bb.map2colormap = function(x, xmax=NA, symmetric=NA, cmap=colorRamp(c('blue','white', 'red')), na_col='cornsilk3'){
  if (is.na(xmax)){
    xmax = max(abs(x), na.rm = T)
  }
  
  x = x/xmax

  # if any value is negative, make symmetric by scaling
  if ((!is.na(symmetric) & symmetric == T) | any(x[!is.na(x)] <0)){
    x = (x+1)/2
  }
  x[x<0 ] =0
  x[x>1] = 1
  
  # apply
  c_x =cmap(x)
  
  # deal with NA
  if (is.character(na_col)){
    na_rgb_col = col2rgb(na_col)
  } else {
    na_rgb_col = na_col
  }
  
  if (any(is.na(c_x))){
    na_fil = which(apply(is.na(c_x),1, any))
    
    for (row in na_fil){
      c_x[row,] = na_rgb_col
    }
  }
  
  rgb_x = rgb(c_x/255)
  return(rgb_x)
}

# calculate percentiles
bb.getPerc <- function(x){
  fkt <- ecdf(x)
  return(fkt(x))
}

bb.censor_dat = function(x, quant = 0.999){
  q = quantile(x,quant)
  x[x>q] = q
  return(x)
  
}

# calculate stats
bb.getStats <- function(df,varName,grpVar,fkt=function(x){x},meltTab = F,bootstrapSD = F){
  df[,tCol := fkt(get(varName))]
  tdt <- df[,list(
    min_c=min(tCol,na.rm=T),
    lower05_c=quantile(tCol, .05, na.rm=TRUE),
    lower25_c=quantile(tCol, .25, na.rm=TRUE),
    median_c=quantile(tCol, .50, na.rm=TRUE),
    mean_c=mean(tCol,na.rm=TRUE),
    mean0.1trim_c=mean(tCol,na.rm=TRUE,trim=0.1),
    upper75_c=quantile(tCol, .75, na.rm=TRUE),
    upper95_c=quantile(tCol, .95, na.rm=TRUE),
    max_c=max(tCol,na.rm=TRUE)),
    by=c(grpVar)]
  
  if (bootstrapSD){
    outBoot = ddply(df,grpVar,function(x){
      
      # mean
      out = boot(x$tCol,statistic= function(x, index) mean(x[index]),R=1000)
      outDat = data.frame(mean_sd_c=sd(out$t))
      allDat = outDat
      
      # median
      out = boot(x$tCol,statistic= function(x, index) median(x[index]),R=1000)
      outDat = data.frame(median_sd_c=sd(out$t))
      allDat = cbind(allDat,outDat)
      
      # trimMn
      out = boot(x$tCol,statistic= function(x, index) mean(x[index],trim=0.1),R=1000)
      outDat = data.frame(mean0.1trim_sd_c=sd(out$t))
      allDat = cbind(allDat,outDat)
      
      return(allDat)
    },.parallel=T)
    
    outBoot = data.table(outBoot)
    setkeyv(outBoot,grpVar)
    setkeyv(tdt,grpVar)
    tdt = tdt[outBoot]
  }
  #delete the temporary column
  df[,tCol:=NULL]
  if (meltTab == T){
    tdt <- melt(tdt,
                id.vars=grpVar,
                variable.factor=F,
                variable.name='stats')
  }
  return(tdt)
}

# plot summary stats
bb.plot_sumStats <- function(df,varName = 'value',
                             condName = 'condition',
                             channelName = 'channel',
                             fkt=function(x){x}){
  stats = bb.getStats(df,varName,c(condName,channelName),fkt,meltTab=T)
  stats = subset(stats,!stats %in% c('max_c','min_c','lower05_c','upper95_c'))
  statsCol = 'stats'
  ycol = 'transfColumn'
  
  p=ggplot(df,aes(x=as.factor(get(condName)),y=fkt(get(varName))), environment = environment())+
    facet_wrap(as.formula(paste("~", channelName)),scale='free',ncol=4)+
    geom_violin()+
    geom_point(data = stats,aes_string(x=paste('as.numeric(as.factor(',condName,'))'),y=varName,colour=statsCol))+
    geom_line(data = stats,aes_string(x=paste('as.numeric(as.factor(',condName,'))'),y=varName,colour=statsCol))
  
  return(p)
}


bb.allComb = function(cDat,idvar='xcat',varvar='ycat',valvar='val'){
  cDat =melt(dcast.data.table(cDat,paste(idvar,varvar,sep='~'),value.var=valvar),
             id.vars = idvar,
             variable.name = varvar,
             value.name = valvar) 
  return(cDat)
}

bb.equalSamp  = function(dat,col='condition'){
  minS = min(dat[, .N,by = get(col)]$N)
  
  dat[,isTaken := sample(c(rep(F,.N-minS),rep(T,minS)),.N),by=get(col)]
  dat = dat[isTaken == T]
  dat[,isTaken :=NULL]
  return(dat)
}



# calculate percentiles
aq.getPerc <- function(x){
  fkt <- ecdf(x)
  return(fkt(x))
}


