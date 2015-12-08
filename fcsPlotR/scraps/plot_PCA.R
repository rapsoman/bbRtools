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


##### Calculate PCA ###

pcadat = prcomp(plotmat[, -'id',with=F], scale. = T)
plot(pcadat)

pca <- as.data.table(pcadat$x)
pca$id = plotmat[, id]
setkey(pca, id)
setkey(dat, id)


ggplot(subset(dat[pca], channel == "Event_length"), aes(x=PC1, y = PC2, col=condition))+
  geom_point(size=4, alpha=1, shape='o')+
  scale_colour_brewer(palette="Paired")

ggplot(subset(dat[pca], channel == "Event_length"), aes(x=PC1, y = PC2, col=bb.censor_dat(counts, 0.99)))+
  geom_point(alpha=0.2, size=2)+
  scale_colour_gradientn(colours=rev(brewer.pal(11,'Spectral')))

# plot as 3d
pdat = subset(dat[pca], channel == "Event_length")

# plot conditions
col = rainbow(5)[as.factor(pdat$condition)]
scatterplot3js(pdat$PC1,pdat$PC3, pdat$PC4, size=0.2, color = col , labels =pdat$condition,stroke=NA )

# plot counts
col = bb.map2colormap(bb.censor_dat(pdat$counts, 0.99), cmap=colorRamp(rev(brewer.pal(11,'Spectral'))))
scatterplot3js(pdat$PC1,pdat$PC2, pdat$PC3, size=0.2, color = col , labels =pdat$condition,stroke=NA )

cor(pdat$PC2, pdat$counts)

### load the pca loadings
pca_loadings = as.data.table(pcadat$rotation)
pca_loadings$channel = rownames(pcadat$rotation)
pca_loadings = melt.data.table(pca_loadings, id.vars = 'channel')

pca_loadings[, pca_dim := as.numeric(gsub("[^0-9]", "",variable))]
ggplot(pca_loadings[ pca_dim < 6], aes(x=channel, y=value))+
  facet_grid(variable~.)+
  geom_bar(stat='identity')

