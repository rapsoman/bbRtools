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
fileDir = '/home/vitoz/imls-bod/Data Stephane/Kinetic Combined'
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
plotHeight =180
# id Columns
idCol = c('Time','Event_length','barcode','condition')

#### Script
condDict = aq.getInfoFromFileList(list.files(fileDir),'_',c(grpVarPos,orderVarPos))
dat <- aq.loadConvertMultiFCS(list.files(fileDir),fileDir,condDict,subSample=subSample)

# remove unecessary columns
dat = dat[,-grep("Beads|BC.|Iridium193|Iridium191|bead|MCB|File\ Number|DNA|Event\ \\#", colnames(dat)), with=FALSE]

# melt
dat =melt(dat,id.vars=idCol,
          variable.factor=F,
          variable.name='channel')

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