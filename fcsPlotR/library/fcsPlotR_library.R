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
aq.loadFCS <-function(filePath){
  datFCS <- read.FCS(filePath,min.limit=NULL,transformation = 'linearize')  
  return(datFCS)
}

aq.flowFrame2dt <-function(datFCS){
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

aq.flowFrame2metaDt <-function(datFCS){
  # gets the metadata table
  metaDt <- data.table(datFCS@parameters@data)
  return(metaDt)
}

aq.getInfoFromString<- function(name,sep='_',strPos=2,censorStr='.fcs'){
  tp <- gsub(censorStr,'',strsplit(name,sep)[[1]][strPos])
  tp <- paste(tp,collapse = sep)
  return(tp)
}

aq.getInfoFromFileList <- function(fileList,sep,strPos,censorStr='.fcs'){
  condDict = sapply(fileList,function(file){aq.getInfoFromString(file,sep,strPos,censorStr)})
  names(condDict) <- fileList
  return(condDict)
}

aq.loadConvertMultiFCS <- function(fileList,fileDir,condDict){
  # fileDir: from which directory should the fcs files be read
  # condition = a named list (dictionary) that assigns filenames to condition
  dat = data.table()
  for(file in fileList){
    if(file_ext(file) =='fcs'){
      cond = condDict[[file]]
      #tp = aq.getInfoFromString(file,sep='_',strPos=condPos,censorStr='.fcs')
      #load and convert to dt
      tmpDT <- aq.flowFrame2dt(aq.loadFCS(file.path(fileDir,file)))
      #combine dt
      dat = rbindlist(list(dat,tmpDT[,'condition':=cond]), fill=T)
    }
  }
  return(dat)
}
### tSNE calculations ####
aq.calcTSNE <- function(input_dat, channels, value_var='counts', channel_var='channel',
                        id_var='id', group_var ='condition', scale=F,
                        subsample_groups=F, subsample_mode='equal',
                        verbose=T,
                        ...){
#Calculates the bhSNE from a melted data table
#input_dat: melted data_table
#channels: list of channels to be Used
#channel_var: column with the channel names
#value_var: column to be used as values
#id_var: column name for the unique cell ID
#scale: should the input data be rescaled? 
#subsample_groups: T: resamples all groups to the group with the least members, Integer: resamples all groups to the integer
#subsample_mode: 
#  'equal': if integer provided to resample_groups all groups will contain min(resample_groups, min_groupsize) cells
#  'unequal': if integer provided to resample_groups, all groups will contain maximally resample_groups cells (or less)
#verbose: should bhSne be verbose?

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
  
  tsne_out <- Rtsne(tsnedat, verbose=verbose,...)
  tsne_out$Y = data.table(bh_1=tsne_out$Y[,1], bh_2=tsne_out$Y[,2], id = dt[, get(id_var)])
  setkeyv(tsne_out$Y, id_var)
  tsne_out$channels = channels
  tsne_out$scale= F
  if (group_var %in% names(input_dat)){
    tsne_out$groups= input_dat[, unique(get(group_var))]
  } else {
    tsne_out$groups = NaN
  }
  return(tsne_out)
  tsne_out$subsample_groups = subsample_groups
  tsne_out$subsample_mode = subsample_mode
}



### various small helper functions ####

# calculate percentiles
aq.getPerc <- function(x){
  fkt <- ecdf(x)
  return(fkt(x))
}

# calculate stats
aq.getStats <- function(df,varName,grpVar,fkt=function(x){x},meltTab = F,bootstrapSD = F){
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
aq.plot_sumStats <- function(df,varName = 'value',
                             condName = 'condition',
                             channelName = 'channel',
                             fkt=function(x){x}){
  stats = aq.getStats(df,varName,c(condName,channelName),fkt,meltTab=T)
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


aq.allComb = function(cDat,idvar='xcat',varvar='ycat',valvar='val'){
  cDat =melt(dcast.data.table(cDat,paste(idvar,varvar,sep='~'),value.var=valvar),
             id.vars = idvar,
             variable.name = varvar,
             value.name = valvar,) 
  return(cDat)
}


