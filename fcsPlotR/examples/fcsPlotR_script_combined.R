library(RColorBrewer)
library(gplots)
library(stringr)

## Load the library
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source('fcsPlotR_library.R')

## Set settings

# folders
#fileDir = '/home/vitoz/imls-bod/Data Stephane/MDM 4x'
#plotDir = '/home/vitoz/imls-bod/Data Stephane/MDM 4x/Plot'
fileDirs = c('/home/vitoz/imls-bod/Data Stephane/kinetic1','/home/vitoz/imls-bod/Data Stephane/kinetic2','/home/vitoz/imls-bod/Data Stephane/kinetic3')
plotDir = '/home/vitoz/imls-bod/Data Stephane/Kinetic Combined/Plots'

# settings
#grpVarPos = 3
#orderVarPos = 1
#doPlot = T
subSample = 10000 # subsampling for a better performance, put NA if not wanted
grpVarPos = 2 # where in the _ separated name is the grouping variable=
orderVarPos = 3 # where in the _ separated name is the variable for sorting within the group?
doPlot = T
plotName = 'Overviewlog10'
fkt = log10
ylabel = 'Conditions'
xlabel = 'log10(Counts)'

# Plot dimensions
plotWidth = 20
plotHeight = 260
# id Columns
idCol = c('Time','Event_length','barcode','condition')

#### Script


dat= data.table()

for (fileDir in fileDirs){
  condDict = aq.getInfoFromFileList(list.files(fileDir),'_',c(grpVarPos,orderVarPos))
  tdat <- aq.loadConvertMultiFCS(list.files(fileDir),fileDir,condDict,subSample=subSample)
  
  # remove unecessary columns
  tdat = tdat[,-grep("Beads|BC.|Iridium193|Iridium191|bead|MCB|File\ Number|DNA|Event\ \\#", colnames(tdat)), with=FALSE]
  
  # melt
  tdat =melt(tdat,id.vars=idCol,
            variable.factor=F,
            variable.name='channel')
  tdat[,channel := paste(channel,substr(fileDir,nchar(fileDir),nchar(fileDir)),sep=' rep ')]
  dat= rbindlist(list(dat,tdat))
}

# censor small values
dat[value < 1,value:=1]

p = aq.plot_sumStats_grouped(dat,varName = 'value',
                     condName = 'condition',
                     channelName = 'channel',
                     fkt=fkt,
                     ylabel=ylabel,
                     xlabel=xlabel)

if (doPlot){
  ggsave(filename = file.path(plotDir,'Overviewlog10_2.pdf'),plot = p,width = plotWidth,height=plotHeight,limitsize=F)
}