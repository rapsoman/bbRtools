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

#### aq.loadConvertMultiFCSaq ####
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

# small helper functions

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




# simulate linear/log-linear/fkt-linear data with real cytof data
# by assuming a maximal mean abundance per cell
# and adding gaussian noise to the predicted data as well as to the
# cytof measurements


aq.simulateData <- function(sumStats,cAssCvSRM=0.1,cAssCvCyt=0.15,cAssMax=1,
                            grpVar='channel',valVar = 'mean_c',fkt=function(x){x},OneToOne=F){
  
  # calc assumed conc and add the gaussian noise to the conc and to the cytof meas
  if(!OneToOne){
    sumStats[ ,
              conc_ass := sapply(cAssMax*fkt(get(valVar))/max(fkt(get(valVar))),
                                 function(x){
                                   rnorm(n = 1,mean =x, sd = x * cAssCvSRM)}),
              by=get(grpVar)]
  } else {
    # If 1to1, the two measurement types are assumed to be related with a direct correspondence.
    sumStats[ ,
              conc_ass := sapply(fkt(get(valVar)),
                                 function(x){
                                   rnorm(n = 1,mean =x, sd = x * cAssCvSRM)}),
              by=get(grpVar)]  
  }
  
  
  sumStats[,mean_ass:= sapply(mean_c,function(x){
    rnorm(n = 1,mean =x, sd = x * cAssCvCyt)})]
  
  return(sumStats)
}


# plot the predictions
aq.predPlot <- function(sumStats,
                        dat,
                        lmForm=y~x,
                        title,
                        xval='mean_ass',
                        yval='conc_ass',
                        condVar = 'condition',
                        pTitle='Prediction accuracy'){
  
  sumStats_perAB = subset(melt(aq.getStats(dat,'value',c('channel')),
                               id.vars=c('channel'),
                               variable.factor=F,
                               variable.name='stats'),!stats %in% c('max_c','min_c') & 
                            channel %in% sumStats$channel)
  
  
  p = ggplot(sumStats,aes_string(x=xval,y=yval))+
    geom_smooth(method='lm',formula=lmForm,fullrange=T)+
    geom_point(aes_string(colour=condVar))+
    scale_colour_discrete(h=c(170,230), 
                          l=rev(seq(0,80,length.out=
                                      sumStats[,length(unique(get(condVar)))]+2))) +
    facet_wrap(~channel,scale='free')+
    expand_limits(x=0,y=0)+
    geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
    #add for legend
    geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
    ggtitle(pTitle)
}

# calculates the predictions with linear models

# plot the predictions
aq.plotMSvsCy<- function(datMS,
                         datCy,
                         lmForm=y~x,
                         ms_val='MS_mn',
                         ms_std = 'MS_std',
                         cy_val='value',
                         cy_summary = 'mean',# which summary stats to use 'mean' or 'median','mean0.1trim'
                         commonID = 'channel',
                         ms_ID ='EG.ProteinId',
                         cy_ID = 'channel',
                         condVar = 'perturbation',
                         groupVar = 'cells',
                         pTitle='MS vs cyTOF',
                         xlab = 'CyTOF [counts/cell]',
                         ylab = 'SRM [copies/cell]',
                         nCols = 4,
                         nPred = 100,
                         sizePoint = 5,
                         fitMethod = 'lm' # lm or weightedLM
){
  
  # merge MS and cyTOF data
  datCy.stats =aq.getStats(datCy,cy_val,unique(c(commonID,cy_ID,condVar,groupVar)),function(x){x},meltTab=F)
  cy_summary2 = paste(cy_summary,'_c',sep='')
  
  setkeyv(datMS,c(commonID,condVar,groupVar))
  setkeyv(datCy.stats,c(commonID,condVar,groupVar))
  
  
  msCols= unique(c(ms_ID,commonID,ms_val,ms_std,condVar,groupVar))
  cyCols = unique(c(cy_ID,commonID,cy_summary2,condVar,groupVar))
  joinDat = merge(datMS[,msCols,with=F],datCy.stats[,cyCols,with=F],all = T,allow.cartesian=TRUE)
  
  joinDat[,groupN := paste(get(cy_ID),'vs',get(ms_ID))]
  joinDat = subset(joinDat,!is.na(get(cy_ID)) & !is.na(get(ms_ID)) & !is.na(get(cy_summary2)) & !is.na(get(ms_val)))
  joinDat[,condition:= paste(cells,perturbation,sep='_')]
  
  
  
  # joinDat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary2),get(ms_val)),sep='\n'),by=groupN]
  #   sumStats_perAB = subset(melt(aq.getStats(dat,'value',c('channel')),
  #                                id.vars=c('channel'),
  #                                variable.factor=F,
  #                                variable.name='stats'),!stats %in% c('max_c','min_c') & 
  #                             channel %in% sumStats$channel)
  
  
  if (fitMethod == 'lm'){
    
    # calculate the equation
    lm_eqn = function(x,y){
      m = lm(y~x);
      eq <- substitute(paste("y=" , a ,"+", b ,"x,\n r2 =",r2), 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      
      return(as.character(eval(eq)));          
    }
    
    joinDat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary2),get(ms_val)),sep='\n'),by=groupN]
    p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
      geom_smooth(method='lm',formula=lmForm,fullrange=T)+
      geom_point(aes_string(colour=groupVar,shape=condVar),size=sizePoint)+
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_wrap(as.formula(paste('~','groupLab')),scale='free',ncol=nCols)+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
  } else if (fitMethod =='weightedLM'){
    
    # fit the  weighted lm
    lmDat = joinDat[,list(
      lmfit= list(lm(get(ms_val)~get(cy_summary2), weights = 1/(get(ms_std)^2))),
      maxX = max(get(cy_summary2)),
      minX = 0
    ),by=groupN]
    
    predDat = lmDat[,list(newX = seq(minX,maxX,length=nPred),
                          lmfit=lmfit),by=groupN]
    setkey(lmDat,groupN)
    setnames(predDat,'newX',cy_summary2)
    predDat = predDat[,list(predY=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD ,interval = "none"),
                            pred_lwr=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD , interval = "confidence", level = 0.95)[,2],
                            pred_upr=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD ,interval = "confidence", level = 0.95)[,3],
                            newX = get(cy_summary2)),
                      by=c('groupN')]
    # 
    setnames(predDat,'newX',cy_summary2)
    setnames(predDat,'predY',ms_val)
    
    #
    p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
      #  geom_smooth(method='lm',formula=lmForm,fullrange=T)+
      geom_point(aes_string(colour=groupVar,shape=condVar),size=sizePoint)+
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_wrap(as.formula(paste('~','groupN')),scale='free',ncol=nCols)+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))+
      geom_line(data = predDat, colour = "blue") + 
      geom_ribbon(mapping = aes(ymax = pred_upr, ymin = pred_lwr), data = predDat, 
                  alpha = 0.4, fill = "grey60")
    
  } else if (fitMethod == 'lmPerCell') {
    
    # calculate the equation
    lm_eqn = function(x,y){
      m = lm(y~x);
      eq <- substitute(paste("y=" , a ,"+", b ,"x,\n r2 =",r2), 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      
      return(as.character(eval(eq)));          
    }
    
    joinDat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary2),get(ms_val)),sep='\n'),by=groupN]
    p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
      
      geom_point(aes_string(colour=groupVar,shape=condVar),size=sizePoint)+
      geom_smooth(method='lm',formula=lmForm,fullrange=F,aes_string(group=groupVar))+
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_wrap(as.formula(paste('~','groupLab')),scale='free',ncol=nCols)+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
    
  }else if (fitMethod == 'lmPerCell2') {
    
    # calculate the equation
    lm_eqn = function(x,y){
      m = lm(y~x);
      eq <- substitute(paste("y=" , a ,"+", b ,"x,\n r2 =",r2), 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      
      return(as.character(eval(eq)));          
    }
    
    joinDat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary2),get(ms_val)),sep='\n'),by=groupN]
    p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
      
      geom_point(aes_string(colour=groupVar,shape=condVar),size=sizePoint)+
      geom_smooth(method='lm',formula=lmForm,fullrange=T,aes_string(group=groupVar))+
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_grid(as.formula(paste('groupLab','~',groupVar)))+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
    
  } else if(fitMethod == 'lmAllPepPlot'){
    
    # calculate the equation
    lm_eqn = function(x,y){
      m = lm(y~x);
      eq <- substitute(paste("y=" , a ,"+", b ,"x,\n r2 =",r2), 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      
      return(as.character(eval(eq)));          
    }
    joinDat[,groupN := paste(get(commonID),'all peptides')]
    joinDat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary2),get(ms_val)),sep='\n'),by=groupN]
    p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
      geom_smooth(method='lm',formula=lmForm,fullrange=T)+
      geom_point(aes_string(colour='condition',shape=ms_ID),size=sizePoint)+
      
      
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_wrap(as.formula(paste('~','groupLab')),scale='free',ncol=nCols)+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
  }
  else  { error('fitMethod not found!')}
  
  
  # setkey(joinDat,groupN)
  # p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
  #   geom_smooth(method='lm',formula='y~x',fullrange=T)+
  #   geom_point(aes_string(colour=groupVar,shape=condVar))+
  # 
  #   #     scale_colour_discrete(h=c(170,230), 
  #   #                           l=rev(seq(0,80,length.out=
  #   #                                       sumStats[,length(unique(get(condVar)))]+2))) +
  #   facet_wrap(as.formula(paste('~',groupVar)),scale='free',ncol=3)+
  #   expand_limits(x=0,y=0)+
  #   # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
  #   #add for legend
  #   #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
  #   ggtitle(pTitle)+
  #   xlab(xlab)+
  #   ylab(ylab)+ 
  #   geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
  # 
  # 
  #   geom_line(data = predDat, colour = "blue") + 
  #   geom_ribbon(mapping = aes(ymax = pred_upr, ymin = pred_lwr), data = predDat, 
  #               alpha = 0.4, fill = "grey60")
  
  
  return(p)
}




aq.fitDualLm = function(joinDat = joinDat, # the data set
                        markCol = 'EG.ModifiedSequence', #
                        condCol ='condition',
                        groupVar = 'EG.ProteinId',
                        lmForm=y~x,
                        xCol='mean_c',
                        yCol='MS_mn',
                        r2Cut=0.2
){
  # This will iteratively adjust the cell counts in order to improve
  # the fit
  
  # fit the lm
  #joinDat[,groupN := paste(get(condCol),'vs',get(markCol))]
  joinDat[,adjustedMS := get(yCol)]
  
  for (i in 1:500){
    setkeyv(joinDat,markCol)
    lmDat = joinDat[,list(
      lmFit= list(lm(adjustedMS~get(xCol)))
    ),by=c(markCol)]
    
    setkeyv(lmDat,markCol)
    # get the r2
    lmDat[,r2 := summary(lmFit[[1]])$r.squared,by=c(markCol)]
    
    # filter by r2
    fitJoinDat=  joinDat[lmDat[r2>r2Cut,get(markCol)]]
    
    #
    
    
    setnames(fitJoinDat,markCol,'markcol')
    # predict values by lm
    fitJoinDat[,predY := predict(lmDat[.BY[1],lmFit][[1]],newdata = .SD),by=markcol]
    fitJoinDat[,condY :=1]
    fitJoinDat[,condX := adjustedMS/predY]
    
    #   # fit the condition specific lm
    #   lmDatCond = fitJoinDat[,list(
    #     lmFit= list(lm(condY ~ condX+0))
    #   ),by=c(condCol)]
    #   
    #   lmDatCond[,corrCoef := coef(lmFit[[1]]),by=c(condCol)]
    #   setkeyv(lmDatCond,condCol)
    #   # predict the adjusted ms values
    #   setkeyv(joinDat,condCol)
    #   joinDat[fitJoinDat[,unique(get(condCol))],adjustedMS := adjustedMS*lmDatCond[get(condCol),corrCoef]]
    setkeyv(joinDat,condCol)
    setkeyv(fitJoinDat,condCol)
    adjDat = fitJoinDat[,list(adjFac = 1/median(condX)),by=c(condCol)]
    setkeyv(adjDat,condCol)
    #joinDat[fitJoinDat[,unique(get(condCol))],adjFac := adjDat[.BY,adjFac],by=c(condCol)]
    
    joinDat[fitJoinDat[,unique(get(condCol))],adjustedMS := adjustedMS*adjDat[.BY,adjFac],by=c(condCol)]
    
  }
  
  joinDat[,adjustedMS := get(yCol)]
  # update each condition randomly by its own
  for (curCond  in sample(rep(joinDat[,unique(get(condCol))],5))){
    setkeyv(joinDat,markCol)
    lmDat = joinDat[,list(
      lmFit= list(lm(adjustedMS~get(xCol)))
    ),by=c(markCol)]
    
    setkeyv(lmDat,markCol)
    # get the r2
    lmDat[,r2 := summary(lmFit[[1]])$r.squared,by=c(markCol)]
    
    # filter by r2
    fitJoinDat=  joinDat[lmDat[r2>r2Cut,get(markCol)]]
    
    #
    
    
    setnames(fitJoinDat,markCol,'markcol')
    # predict values by lm
    fitJoinDat[,predY := predict(lmDat[.BY[1],lmFit][[1]],newdata = .SD),by=markcol]
    fitJoinDat[,condY :=1]
    fitJoinDat[,condX := adjustedMS/predY]
    
    #   # fit the condition specific lm
    #   lmDatCond = fitJoinDat[,list(
    #     lmFit= list(lm(condY ~ condX+0))
    #   ),by=c(condCol)]
    #   
    #   lmDatCond[,corrCoef := coef(lmFit[[1]]),by=c(condCol)]
    #   setkeyv(lmDatCond,condCol)
    #   # predict the adjusted ms values
    #   setkeyv(joinDat,condCol)
    #   joinDat[fitJoinDat[,unique(get(condCol))],adjustedMS := adjustedMS*lmDatCond[get(condCol),corrCoef]]
    setkeyv(joinDat,condCol)
    setkeyv(fitJoinDat,condCol)
    adjDat = fitJoinDat[,list(adjFac = 1/mean(condX,trim=0.2)),by=c(condCol)]
    setkeyv(adjDat,condCol)
    #joinDat[fitJoinDat[,unique(get(condCol))],adjFac := adjDat[.BY,adjFac],by=c(condCol)]
    
    joinDat[curCond,adjustedMS := adjustedMS*adjDat[.BY,adjFac],by=c(condCol)]
    
  }
  nCols =7
  pTitle = 'bla'
  p2 = ggplot(joinDat,aes_string(x=xCol,y='adjustedMS'))+
    geom_smooth(method='lm',formula=lmForm,fullrange=T)+
    geom_point(aes_string(colour=condCol))+
    #     scale_colour_discrete(h=c(170,230), 
    #                           l=rev(seq(0,80,length.out=
    #                                       sumStats[,length(unique(get(condVar)))]+2))) +
    geom_point(data=joinDat,aes_string(x=xCol,y=yCol))+
    facet_wrap(as.formula(paste('~',markCol)),scale='free',ncol = nCols)+
    expand_limits(x=0,y=0)
  # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
  #add for legend
  #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
  #  ggtitle(pTitle)+
  #  xlab(xlab)+
  #  ylab(ylab)+
  #  geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
  return(p2)
  
  
  #   setkey(lmDat,groupN)
  #   setnames(predDat,'newX',cy_summary2)
  #   predDat = predDat[,list(predY=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD ,interval = "none"),
  #                           pred_lwr=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD , interval = "confidence", level = 0.95)[,2],
  #                           pred_upr=predict(lmDat[.BY]$lmfit[[1]],newdata = .SD ,interval = "confidence", level = 0.95)[,3],
  #                           newX = get(cy_summary2)),
  #                     by=c('groupN')]
  #   
  #   setnames(predDat,'newX',cy_summary2)
  #   setnames(predDat,'predY',ms_val)
  
  
  
}


aq.allComb = function(cDat,idvar='xcat',varvar='ycat',valvar='val'){
  cDat =melt(dcast.data.table(cDat,paste(idvar,varvar,sep='~'),value.var=valvar),
             id.vars = idvar,
             variable.name = varvar,
             value.name = valvar,) 
  return(cDat)
}


# plot the predictions
aq.plotMSvsCy_complete<- function(dat,
                                  lmForm=y~x,
                                  ms_val='ms.copynr.per.cell',
                                  ms_summary='MS_mn',
                                  ms_summary_sd = 'MS_std_comb',
                                  
                                  ms_ID = 'EG.StrippedSequence',
                                  cy_ID = 'channel',
                                  cy_summary = "cytof_mean0.1trim_c",# which summary stats to use 'mean' or 'median','mean0.1trim'
                                  cy_summary_sd = 'cytof_mean0.1trim_sd_c',
                                  condVar = 'perturbation',
                                  groupVar = 'cells',
                                  pTitle='MS vs cyTOF',
                                  xlab = 'CyTOF [counts/cell]',
                                  ylab = 'SRM [copies/cell]',
                                  nCols = 4,
                                  nPred = 100,
                                  sizePoint = 5,
                                  fitMethod = 'lm' # lm or weightedLM
){
  
  # merge MS and cyTOF data
  
  
  if (fitMethod == 'lm'){
    
    # calculate the equation
    lm_eqn = function(x,y){
      m = lm(y~x);
      eq <- substitute(paste("y=" , a ,"+", b ,"x,\n r2 =",r2), 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      
      return(as.character(eval(eq)));          
    }
    dat[,groupN := paste(get(cy_ID),'vs',get(ms_ID))]
    dat[,groupLab  := paste(groupN,lm_eqn(get(cy_summary),get(ms_summary)),sep='\n'),by=groupN]
    p = ggplot(dat,aes_string(x=cy_summary,y=ms_summary))+
      geom_smooth(method='lm',formula=lmForm,fullrange=T)+
      geom_point(aes_string(x=cy_summary,y=ms_val,colour=groupVar,shape=condVar),size=sizePoint)+
      #     scale_colour_discrete(h=c(170,230), 
      #                           l=rev(seq(0,80,length.out=
      #                                       sumStats[,length(unique(get(condVar)))]+2))) +
      facet_wrap(as.formula(paste('~','groupLab')),scale='free',ncol=nCols)+
      expand_limits(x=0,y=0)+
      # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
      #add for legend
      #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
      ggtitle(pTitle)+
      xlab(xlab)+
      ylab(ylab)+
      geom_errorbar(aes_string(ymax=paste(ms_summary,'+ ',ms_summary_sd),ymin=paste(ms_summary,'-',ms_summary_sd),
                               xmax=paste(cy_summary,'+ ',cy_summary_sd),xmin=paste(cy_summary,'-',cy_summary_sd)))+
      geom_errorbarh(aes_string(xmax=paste(cy_summary,'+ ',cy_summary_sd),xmin=paste(cy_summary,'-',cy_summary_sd)))
  }   else  { error('fitMethod not found!')}
  
  
  # setkey(joinDat,groupN)
  # p = ggplot(joinDat,aes_string(x=cy_summary2,y=ms_val))+
  #   geom_smooth(method='lm',formula='y~x',fullrange=T)+
  #   geom_point(aes_string(colour=groupVar,shape=condVar))+
  # 
  #   #     scale_colour_discrete(h=c(170,230), 
  #   #                           l=rev(seq(0,80,length.out=
  #   #                                       sumStats[,length(unique(get(condVar)))]+2))) +
  #   facet_wrap(as.formula(paste('~',groupVar)),scale='free',ncol=3)+
  #   expand_limits(x=0,y=0)+
  #   # geom_vline(data=sumStats_perAB,aes(xintercept=value,linetype=stats),alpha=0.8)+
  #   #add for legend
  #   #geom_line(data=sumStats_perAB,aes(x=value,y=Inf,linetype=stats))+
  #   ggtitle(pTitle)+
  #   xlab(xlab)+
  #   ylab(ylab)+ 
  #   geom_errorbar(aes_string(ymax=paste(ms_val,'+ ',ms_std),ymin=paste(ms_val,'-',ms_std)))
  # 
  # 
  #   geom_line(data = predDat, colour = "blue") + 
  #   geom_ribbon(mapping = aes(ymax = pred_upr, ymin = pred_lwr), data = predDat, 
  #               alpha = 0.4, fill = "grey60")
  
  
  return(p)
}

