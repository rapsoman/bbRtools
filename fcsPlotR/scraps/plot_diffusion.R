source('/mnt/imls-bod/LabCode/bbRtools/fcsPlotR/library/fcsPlotR_library.R')
library('data.table')
library('RColorBrewer')
library('destiny')
library('plot3D')
library(threejs)
#### Plot fcs files in a folder as a bhSNE ####
# Copy this script in your folder
#


# input files
fcs_folder = '/mnt/imls-bod/LabCode/bbRtools/ExampleData/FCStest/'

# define exlcluded crap channels
crap_channels = c("Time", "Event_length", "MCB102", "MCB104", "MCB105", "MCB106", "MCB108","MCB110","MCB113","MCB115" ,"Beads140" ,
                  "DNA193",    "DNA191" ,   "Live194",      "Live195" ,     "beadDist",     "barcode"
)


########## Start script ###
# load data
fcs_files = list.files(fcs_folder)
dat = bb.loadConvertMultiFCS(fileDir = fcs_folder)

# give each cell an own ID
dat[, id := paste(1:.N, .BY), by=condition]
setkey(dat, id)
# 'melt' data: make a column 'channel' with all the channels and a column 'counts' with all the counts
dat = melt.data.table(dat, id.vars=c('condition','id'), variable.name='channel', value.name = 'counts' , variable.factor = FALSE)


### calculae transformed counts ####
dat[ , counts_transf := asinh(counts/5)]


# make a list
good_channels = unique(dat$channel)[!unique(dat$channel) %in% crap_channels]
print('Channels used:')
print(good_channels)

plotmat = dcast.data.table(subset(dat, channel %in% good_channels), formula = id ~ channel, value.var = 'counts_transf')
###### Run diffusion maps #####



plotmat = plotmat[ id %in% sample(id, 25000),]

difout <- DiffusionMap(scale(plotmat[, -'id', with=F]), verbose = T)

difdat = as.data.table(as.data.frame(difout))
# clean out the unneeded data
difdat <- difdat[, -colnames(plotmat[, -'id', with=F]), with=F]
# set id
difdat[, id:= plotmat$id]
setkey(difdat, id)


                                                                                                                                                                                          ggplot(subset(dat[difdat], channel == "Slamf7"), aes(x=PC1, y = PC2, col=condition))+
  geom_point(size=4, alpha=1, shape='o')+
  scale_colour_brewer(palette="Paired")

ggplot(subset(dat[difdat], channel == "Slamf7"), aes(x=DC1, y = DC2, col=bb.censor_dat(counts, 0.99)))+
  geom_point()+
  scale_colour_gradientn(colours=rev(brewer.pal(11,'Spectral')))

# plot as 3d
pdat = subset(dat[difdat], channel == "Slamf7")

# plot conditions
col = rainbow(5)[as.factor(pdat$condition)]
scatterplot3js(pdat$PC1,pdat$DC1, pdat$DC2, size=0.2, color = col , labels =pdat$condition,stroke=NA )

# plot counts
col = bb.map2colormap(bb.censor_dat(pdat$counts, 0.99), cmap=colorRamp(rev(brewer.pal(11,'Spectral'))))
scatterplot3js(pdat$DC1,pdat$DC2, pdat$DC3, size=0.2, color = col , labels =pdat$condition,stroke=NA )



