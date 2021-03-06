#### Helper functions for loading FCS files ####


#'  Wrapper for flow Core read.FCS
#'
#' @param filePath path to the FCS file
#' @return A flowframe object
#' @importFrom flowCore read.FCS
#' @export
loadFCS <-function(filePath, ...){
  # a wrapper to the standard read.FCS, optimized for cyTOF data
  datFCS <- flowCore::read.FCS(filePath,min.limit=NULL,transformation = FALSE, truncate_max_range = FALSE, ...)
  return(datFCS)
}

#' Converts FlowFrame to data table
#'
#' @param datFCS A flowframe e.g. from read.FCS or loadFCS
#' @return A data.table with row=cell, column=Channel
#' @import data.table
#' @importFrom flowCore exprs
flowFrame2dt <-function(datFCS){
  # converts a flow frame from read.FCS to a data table
  dt = flowCore::exprs(datFCS)
  colnames(dt) = make.names(make.unique(colnames(dt)))
  fil = !is.na(datFCS@parameters@data$desc)
  
  if (any(fil)){
    colnames(dt[,fil]) <- datFCS@parameters@data$desc[fil]
  }
  
  dt <- data.table::data.table(dt)
  if (any(fil)){
    uni_newnames = make.names(make.unique(datFCS@parameters@data$desc))
    data.table::setnames(dt,colnames(dt)[fil],uni_newnames[fil])
  }
  
  colnames(dt) <- make.names(make.unique(colnames(dt)))
  return(dt)
}

#' Loads multiple FCS files at once
#'
#' @param fileList Vector of filename to be loaded. If NaN all fcs files are loaded
#' @param fileDir the base directory for the files or a directory with FCS files
#' @param condDict a named vector which maps Filenames to condition
#'
#' @return a data/tab;e with all combined FCS files
#' 
#' This is a good way to load a whole directory of files (if fileList=NaN).
#' If a condition dictionary (condDict) is provided, the 'condition' column will
#' be equal to the conditions. Otherwise it will be the filename used.
#' @export
#' @import data.table
loadConvertMultiFCS <- function(fileList=NaN,fileDir=NaN,
                                condDict=NaN, subSample = NA,
                                subSampleMode = 'number',
                                preprocessfkt=NULL){
  if (is.na(fileList)){
    fileList = list.files(fileDir, pattern='*.fcs')
  }
  dat = data.table::data.table()
  for(file in fileList){
    if(tools::file_ext(file) =='fcs'){
      if (!is.na(condDict)){
        cond = condDict[[file]]
      } else {
        cond = file
      }
      
      #load and convert to dt
      tmpDT <- flowFrame2dt(loadFCS(file.path(fileDir,file)))
      
      if (!is.na(subSample) && is.numeric(subSample)){
        tmpDT = tmpDT[subsample_unequal(1:nrow(tmpDT), n=subSample, mode=subSampleMode), ]
      }
      if (!is.null(preprocessfkt)){
        tmpDT <- preprocessfkt(tmpDT)
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
#flowFrame2metaDt <-function(datFCS){
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
#' @export
getInfoFromString<- function(name,sep='_',strPos=2,censorStr='.fcs'){
  tp <- gsub(censorStr,'',strsplit(name,sep)[[1]][strPos])
  tp <- paste(tp,collapse = sep)
  return(tp)
}


#' Parses a list of strings and gets the corresponding information
#'
#' See getInfoFromString for more information
#'
#' @export
getInfoFromFileList <- function(fileList,sep='_',strPos=2,censorStr='.fcs'){
  condDict = sapply(fileList,function(file){getInfoFromString(file,sep,strPos,censorStr)})
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
#'  'fraction': subsamples a fraction of the cells to get the total amount of cells indicated in subsample_groups
#' @param verbose should bhSne be verbose?
#' @param dimsd number of dimensions
#' @return tsne object
#' @export
#' @import data.table
#' @import Rtsne.multicore
#' @import Rtsne
calcTSNE <- function(input_dat, channels, value_var='counts', channel_var='channel',
                     id_var='id', group_var ='condition', scale=F,
                     subsample_groups=F, subsample_mode='equal',
                     verbose=T,
                     dims=2, multicore=NULL,
                     ...){
  
  # do subsampling
  if (is.numeric(subsample_groups)){
    if (subsample_mode == 'equal'){
      min_n = min(input_dat[get(channel_var) == channels[1], list(n= .N), by=get(group_var)]$n)
      subsample_groups = min(subsample_groups, min_n)
      ids = input_dat[get(channel_var) == channels[1], list(fil = get(id_var)[sample.int(.N, min(subsample_groups, .N))]),by=get(group_var)]$fil
    } else if (subsample_mode == 'unequal'){
      ids = input_dat[get(channel_var) == channels[1], list(fil = get(id_var)[sample.int(.N, min(subsample_groups, .N))]),by=get(group_var)]$fil
    } else if (subsample_mode == 'fraction'){
      frac_n= input_dat[!duplicated(get(id_var)), subsample_groups/.N]
      if (frac_n > 1){
        warning('no subsampling required as less than requested cells available')
        frac_n = 1
      }
      ids = input_dat[!duplicated(get(id_var)), list(id = sample(get(id_var), ceiling(.N*frac_n))), by = get(group_var)]$id
    } else {
      stop('subsample_mode not a valid option!')
    }
    
  } else if (subsample_groups){
    min_n = min(input_dat[get(channel_var) == channels[1],list(n= .N), by=get(group_var)]$n)
    ids = input_dat[get(channel_var) == channels[1], list(fil = get(id_var)[sample.int(.N, min_n)]),by=get(group_var)]$fil
    
  } else {
    ids = input_dat[!duplicated(get(id_var)), get(id_var)]
  }
  
  
  
  dt = data.table::dcast(input_dat[(get(channel_var) %in% channels) & get(id_var) %in% ids], formula = stats::as.formula(paste(id_var, '~', channel_var)), value.var = value_var)
  
  if (scale){
    tsnedat = scale(dt[, channels, with=F])
    
  } else {
    tsnedat = dt[, channels, with=F]
  }
  
  if (is.null(multicore)){
    tsne_out <- Rtsne::Rtsne(tsnedat, verbose=verbose,dims=dims,...)
  } else {
    tsne_out <- Rtsne.multicore::Rtsne.multicore(tsnedat, verbose=verbose,dims=dims, num_threads=multicore, ...)
  }
  tsne_out$Y = data.table(tsne_out$Y)
  setnames(tsne_out$Y, names(tsne_out$Y), paste("bh",names(tsne_out$Y),sep='_'))
  tsne_out$Y$id = dt[, get(id_var)]
  data.table::setnames(tsne_out$Y, 'id', id_var)
  data.table::setkeyv(tsne_out$Y, id_var)
  tsne_out$channels = channels
  tsne_out$scale= F
  tsne_out$multicore = multicore
  if (group_var %in% names(input_dat)){
    tsne_out$groups= input_dat[, unique(get(group_var))]
  } else {
    tsne_out$groups = NaN
  }
  
  tsne_out$subsample_groups = subsample_groups
  tsne_out$subsample_mode = subsample_mode
  return(tsne_out)
  
}

subsample_unequal <- function(ids, n, mode='number'){
  nids = length(ids)
  if (mode == 'fraction'){
    n = floor(nids * n)
  }
  n = min(n, nids)
  sampids = sample(ids, n, replace=F)
  return(sampids)
}



### Network plotting functions ####

#' Loads a Cytoscape JSON as an igraph graph
#'
#' @param  jdat a loaded Cytoscape JSON
#' @return an igraph-graph
#'
#' @examples
#' 
#' \dontrun{
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
#'}
#' @export
#' @importFrom igraph graph_from_edgelist set.vertex.attribute set.edge.attribute V
load_json_graph = function(jdat){
  "
  Loads a graph with metadata from a cytoscape.js JSON file
  "
  edgelist = sapply(jdat$elements$edges, function(e){
    # Loads the edge matrix f
    return(c(e$data$source, e$data$target))
  })
  # Create graph
  g = igraph::graph_from_edgelist(t(edgelist))
  
  # Loop over the nodes to get metadata
  for (v in jdat$elements$nodes){
    
    idx = which(igraph::V(g)$name == v$data$id)
    
    # Loop over the metadata entries
    for (at in names(v$data)){
      
      g = igraph::set.vertex.attribute(g, at, idx, v$data[[at]])
    }
    
    # set position metadata
    g = igraph::set.vertex.attribute(g,'posx',idx, as.numeric(v$position$x))
    g = igraph::set.vertex.attribute(g,'posy',idx, as.numeric(v$position$y))
    
  }
  
  # Loop over edges metadata
  for (idx in 1:length(jdat$elements$edges)){
    e = jdat$elements$edges[[idx]]
    for (at in names(e$data)){
      g = igraph::set.edge.attribute(g,at,idx, e$data[[at]])
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
#' @export
map2colormap = function(x, xmax=NA, symmetric=NA, cmap=colorRamp(c('blue','white', 'red')), na_col='cornsilk3'){
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

#' calculate percentiles
#' For each value in x  calculate which percentile it is 
#'
#' @param x values to be mapped
#'
#' @return returns the percentile of each value of x
#' @export
getPerc <- function(x){
  fkt <- ecdf(x)
  return(fkt(x))
}

#' censor dat
#' remove the outliers on the upper side by capping the values at the provided quantile
#'
#' @param x values to censor
#' @param quant quantile to censor, i.e. how many percent of values are considered outliers
#' @param symmetric censor on both side. In this case the outliers are assumed to be symetric on both sides. For example if a quantile of 5\% (0.05) is choosen, in the symetric case 2.5\% (0.025) of values are censored on both sides.
#'
#' @return returns the percentile of each value of x
#' @export
censor_dat = function(x, quant = 0.999, symmetric=F){
  if (symmetric){
    lower_quant = (1-quant)/2
    quant = quant+lower_quant
  }
  q = stats::quantile(x,quant)
  x[x>q] = q
  
  if(symmetric){
    q = stats::quantile(x, lower_quant)
    x[x < q] = q
  }
  return(x)
}

#' calculate stats
#' @export
#' @import data.table
getStats <- function(df,varName,grpVar,fkt=function(x){x},meltTab = F,bootstrapSD = F){
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
    outBoot = plyr::ddply(df,grpVar,function(x){
      
      # mean
      out = boot::boot(x$tCol,statistic= function(x, index) mean(x[index]),R=1000)
      outDat = data.frame(mean_sd_c=sd(out$t))
      allDat = outDat
      
      # median
      out = boot::boot(x$tCol,statistic= function(x, index) median(x[index]),R=1000)
      outDat = data.frame(median_sd_c=sd(out$t))
      allDat = cbind(allDat,outDat)
      
      # trimMn
      out = boot::boot(x$tCol,statistic= function(x, index) mean(x[index],trim=0.1),R=1000)
      outDat = data.frame(mean0.1trim_sd_c=sd(out$t))
      allDat = cbind(allDat,outDat)
      
      return(allDat)
    },.parallel=T)
    
    outBoot = data.table::data.table(outBoot)
    data.table:setkeyv(outBoot,grpVar)
    data.table:setkeyv(tdt,grpVar)
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

#' plot summary stats
#' @export
#' @import ggplot2
plot_sumStats <- function(df,varName = 'value',
                          condName = 'condition',
                          channelName = 'channel',
                          fkt=function(x){x}){
  stats = getStats(df,varName,c(condName,channelName),fkt,meltTab=T)
  stats = subset(stats,!stats %in% c('max_c','min_c','lower05_c','upper95_c'))
  statsCol = 'stats'
  ycol = 'transfColumn'
  
  p=ggplot(df,aes(x=as.factor(get(condName)),y=fkt(get(varName))), environment = environment())+
    facet_wrap(as.formula(paste("~", channelName)),scales = 'free', ncol=4)+
    geom_violin()+
    geom_point(data = stats,aes_string(x=paste('as.numeric(as.factor(',condName,'))'),y=varName,colour=statsCol))+
    geom_line(data = stats,aes_string(x=paste('as.numeric(as.factor(',condName,'))'),y=varName,colour=statsCol))
  
  return(p)
}

#' @export
allComb = function(cDat,idvar='xcat',varvar='ycat',valvar='val'){
  cDat = data.table::melt(data.table::dcast(cDat,paste(idvar,varvar,sep='~'),value.var=valvar),
                          id.vars = idvar,
                          variable.name = varvar,
                          value.name = valvar)
  return(cDat)
}

#' @export
#' @import data.table
equalSamp  = function(dat,col='condition'){
  minS = min(dat[, .N,by = get(col)]$N)
  
  dat[,isTaken := sample(c(rep(F,.N-minS),rep(T,minS)),.N),by=get(col)]
  dat = dat[isTaken == T]
  dat[,isTaken :=NULL]
  return(dat)
}


#' Calculate a correlation matrix from a data table
#' @export get_cormat
#' @import data.table
#' @import dplyr
get_cormat <- function(data, xcol, ycol, valuecol, method='pearson', pval = F){
  #' @param data: a data table
  #' @param xcol: the column name to be taken as x column
  #' @param ycol: the column name to be taken as y column
  #' @param valuecol: the column name to be taken as valuecol
  #' @param method: pearson or spearman, type argumetn from 'rcorr'
  #' @param pval: return p values for the correlations?
  
  cormat = dcast.data.table(data, paste(xcol, ycol,sep='~'), value.var = valuecol) %>%
    dplyr::select(-matches(xcol)) %>%
    as.matrix() %>%
    rcorr(type=method)
  if (pval == F){
    cormat = cormat$r
  }
  return(cormat)
}

#'  Run flowsom on a melted data table
#' @export do_flowsom
#' @import data.table
#' @importFrom flowCore flowFrame
#' @import ConsensusClusterPlus
#' @import FlowSOM
do_flowsom <- function(data, channels, valuevar= 'counts_transf', 
                       channelvar='channel', idvar='id', k=20, seed=FALSE,
                       subsample=FALSE, return_output=FALSE, ...){
  #' @param data a data frame in the long format
  #' @param channels a list of channel names to use
  #' @param valuevar the column to take as value variable
  #' @param channelvar the column to take as channel variable
  #' @param idvar the column to toake as id variable
  #' @param k the number of clusters to use for metaclustering
  #' @param seed the random seed - if provided
  #' @param subsample number of cells to use for subsampling - if provided
  #' @param return_output should the cluster algorithm output be provided? If so the return value will be a list with the table and the original output.
  #' @return at table with a column cluster and a column idvar
  
  
  pheno_dat = data.table::dcast.data.table(data[get(channelvar) %in% channels, ],paste(idvar, channelvar,sep='~'),value.var = valuevar)
  all_ids = pheno_dat[, get(idvar)]
  
  if (subsample == FALSE){
    subsample=1
  } 
  
  sampids = pheno_dat[, sample(get(idvar), floor(.N*subsample),replace = F)]
  pheno_dat_samp = pheno_dat[get(idvar) %in% sampids, ]
  ids = pheno_dat_samp[, get(idvar)]
  pheno_dat_samp[, (idvar):=NULL]
  set.seed(seed)
  
  data_FlowSOM <- flowCore::flowFrame(as.matrix(pheno_dat_samp))
  
  
  ###################
  ### RUN FLOWSOM ###
  ###################
  
  # set seed for reproducibility
  
  if (seed) set.seed(seed)
  
  # run FlowSOM (initial steps prior to meta-clustering)
  
  out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out, ...)
  out <- FlowSOM::BuildMST(out)
  
  # extract cluster labels (pre meta-clustering) from output object
  
  labels_pre <- out$map$mapping[, 1]
  
  
  # run meta-clustering
  
  # note: In the current version of FlowSOM, the meta-clustering function 
  # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
  # correctly, so results are not reproducible. We use the internal function 
  # ConsensusClusterPlus::ConsensusClusterPlus() to get around this. However, this
  # will be fixed in the next update of FlowSOM (version 1.5); then the following 
  # (simpler) code can be used instead:
  
  if (seed) {
    out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed,writeTable = F,plot = 'pngBMP')
    out <- out[[k]]$consensusClass
  } else {
    out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k)
  }
  
  # extract cluster labels from output object
  
  cluster <- out[labels_pre]
  
  pheno_clust = data.table::data.table(cluster)
  pheno_clust[, (idvar):=ids]
  pheno_clust[, cluster:=factor(cluster)]
  data.table::setkeyv(pheno_clust, idvar)
  
  if (return_output== FALSE){
    return(pheno_clust)
  } else {
    return(list(pheno_clust, out))
  }
}

#' Make a phenograph from a melted dataframe
#' @export do_phenograph
#' @import data.table
#' @importFrom Rphenograph Rphenograph
#' @importFrom igraph membership
do_phenograph<- function(data, channels, valuevar= 'counts_transf', channelvar='channel',
                         idvar='id', k=20, seed=FALSE, subsample=FALSE, return_output=FALSE, ...){
  #' @param data a data frame in the long format
  #' @param channels a list of channel names to use
  #' @param valuevar the column to take as value variable
  #' @param channelvar the column to take as channel variable
  #' @param idvar the column to toake as id variable
  #' @param k the number of nearest neighbours
  #' @param seed the random seed - if provided
  #' @param subsample number of cells to use for subsampling - if provided
  #' @param return_output should the cluster algorithm output be provided? If so the return value will be a list with the table and the original output.
  #' @return at table with a column cluster and a column idvar
  
  
  
  pheno_dat = data.table::dcast.data.table(data[get(channelvar) %in% channels, ],paste(idvar, channelvar,sep='~'),value.var = valuevar)
  
  all_ids = pheno_dat[, get(idvar)]
  
  if (subsample == FALSE){
    subsample=1
  } 
  if (seed){
    set.seed(seed)
  }
  sampids = pheno_dat[, sample(get(idvar), floor(.N*subsample),replace = F)]
  pheno_dat_samp = pheno_dat[get(idvar) %in% sampids, ]
  ids = pheno_dat_samp[, get(idvar)]
  pheno_dat_samp[, (idvar):=NULL]
  
  rpheno_out = cytofkit::Rphenograph(pheno_dat_samp, k, ...)
  cluster = igraph::membership(rpheno_out)
  names(cluster) <- as.character(as.integer(names(cluster)))
  pheno_clust = data.table::data.table(x=ids)
  setnames(pheno_clust, 'x', idvar)
  pheno_clust[, cluster:=factor(cluster[as.character(seq(length(ids)))])]
  data.table::setkeyv(pheno_clust, idvar)
  
  if (return_output == FALSE){
    return(pheno_clust)
  } else {
    # Rphenograph has changed the output to return [igraph_graph, communities]
    # Before it was only communities. To not break compatiblities, the order of this
    # output is changed.
    return(list(pheno_clust, rpheno_out[[2]], rpheno_out[[1]]))
  }
}
