#### Plot fcs files in a folder as a bhSNE ####

# initialize
# 
source('fcsPlotR/library/fcsPlotR_library.R')
library('data.table')
library('Rtsne')

# input files
fcs_folder = '/mnt/imls-bod/for Stephane/FCS test/'

# load data
fcs_files = list.files(fcs_folder)
cond_dict = fcs_files
names(cond_dict) = fcs_files
dat = aq.loadConvertMultiFCS(fileList = fcs_files, fileDir = fcs_folder, condDict = cond_dict)

# give each cell an own ID
dat[, id := 1:.N]
setkey(dat, id)
# 'melt' data: make a column 'channel' with all the channels and a column 'counts' with all the counts
dat = melt.data.table(dat, id.vars=c('condition','id'), variable.name='channel', value.name = 'counts' , variable.factor = FALSE)


### calculae transformed counts ####
dat[ , counts_transf := asinh(counts/5)]

###### Run Tsne #####

# define exlcluded crap channels
crap_channels = c("Time", "Event_length", "MCB102", "MCB104", "MCB105", "MCB106", "MCB108","MCB110","MCB113","MCB115" ,"Beads140" ,
                  "DNA193",    "DNA191" ,   "Live194",      "Live195" ,     "beadDist",     "barcode"
                  )
# make a list
good_channels = unique(dat$channel)[!unique(dat$channel) %in% crap_channels]
print(good_channels)

# make a matrix with all channels used for tsne
# use transformed counts

mat =as.matrix(dcast.data.table(subset(dat, channel %in% good_channels), formula = id ~ channel, value.var = 'counts_transf'))
ids = mat[,1]
mat = mat[,-1]


# a) run without scaling 
tsne_out = Rtsne(mat)
bh_dat = data.table(bh_1_unscaled=tsne_out$Y[,1], bh_2_unscaled=tsne_out$Y[,2], id = ids)

setkey(bh_dat, id)
setkey(dat, id)
dat = dat[bh_dat]

# b) run with scaling
tsne_out = Rtsne(scale(mat))

tsne_out$Y = data.table(bh_1_scaled=tsne_out$Y[,1], bh_2_scaled=tsne_out$Y[,2], id = ids)

setkey(bh_dat, id)
setkey(dat, id)
dat = dat[bh_dat]



p = ggplot(subset(dat, !duplicated(id)), aes(x=bh_1_unscaled, y=bh_2_unscaled))+
  geom_point(alpha=0.2, size=0.5)
p

p = ggplot(subset(dat, !duplicated(id)), aes(x=bh_1_scaled, y=bh_2_scaled, color=counts))+
  geom_point(alpha=0.3, size=0.5)
p

good_channels
marker = good_channels[1]

censor_dat = function(x, quant = 0.999){
  q = quantile(x,quant)
  x[x>q] = q
  return(x)
  
}

p = ggplot(subset(dat, channel == marker), aes(x=bh_1_scaled, y=bh_2_scaled, color=censor_dat(counts)))+
  geom_point(alpha=1, size=1)+
  ggtitle(marker)
p

pdf(file.path(fcs_folder,'test_1.pdf'),width = 12, height = 8)

p = ggplot(subset(dat, !duplicated(id)), aes(x=bh_1_scaled, y=bh_2_scaled, color=condition))+
  geom_point(alpha=1, size=0.5)+
  scale_colour_brewer(palette="Paired")+
  ggtitle('Conditions')+
  guides(color=guide_legend(override.aes=list(size=5)))
print(p)

for (chan in unique(dat$channel)){
  p = ggplot(subset(dat, channel == chan), aes(x=bh_1_scaled, y=bh_2_scaled, color=censor_dat(counts)))+
    geom_point(alpha=1, size=0.5)+
    ggtitle(chan)
  print(p)
}



dev.off()
