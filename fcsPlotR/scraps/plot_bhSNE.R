source('/mnt/imls-bod/LabCode/bbRtools/fcsPlotR/library/fcsPlotR_library.R')
library('data.table')
library('Rtsne')
library(threejs)
library('RColorBrewer')
library(plotly)
#### Plot fcs files in a folder as a bhSNE ####
# Copy this script in your folder
#
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

###### Run Tsne #####


# make a list
good_channels = unique(dat$channel)[!unique(dat$channel) %in% crap_channels]
print(good_channels)
setkey(dat, id)

# make a matrix with all channels used for tsne
tsne_out_unscaled = bb.calcTSNE(dat, good_channels, value_var = 'counts_transf', 
                                id_var = 'id', group_var='condition', scale = F, subsample_groups=F, subsample_mode='equal')
tsne_out_scaled = bb.calcTSNE(dat, good_channels, value_var = 'counts_transf', id_var = 'id', scale = T,
                               group_var='condition',  subsample_groups=F, subsample_mode='equal')

tsne_out_scaled_3d = bb.calcTSNE(dat, good_channels, value_var = 'counts_transf', id_var = 'id', scale = T,
                              group_var='condition',  subsample_groups=F, subsample_mode='equal',  dims=3)

tsne_out_unscaled_sub = bb.calcTSNE(dat, good_channels, value_var = 'counts_transf', 
                                id_var = 'id', group_var='condition', scale = F, subsample_groups=T, subsample_mode='equal')
tsne_out_scaled_sub = bb.calcTSNE(dat, good_channels, value_var = 'counts_transf', id_var = 'id', scale = T,
                              group_var='condition',  subsample_groups=T, subsample_mode='equal')


p = ggplot(subset(dat, channel == "Event_length")[tsne_out_scaled$Y], aes(x=bh_V1, y=bh_V2, color=condition))+
  geom_point(alpha=0.2, size=2)
p

p = ggplot(subset(dat, channel == "Event_length")[tsne_out_scaled$Y], aes(x=bh_V1, y=bh_V2, color=bb.censor_dat(counts, 0.99)))+
  geom_point(alpha=0.2, size=2)+
  scale_colour_gradientn(colours=rev(brewer.pal(11,'Spectral')))
p


p = ggplot(subset(dat, channel == "Event_length")[tsne_out_unscaled$Y], aes(x=bh_V1, y=bh_V2, color=condition))+
  geom_point(alpha=0.3, size=2)
p

p = ggplot(subset(dat, channel == "Event_length")[tsne_out_unscaled_sub$Y], aes(x=bh_V1, y=bh_V2, color=condition))+
  geom_point(alpha=0.3, size=2)
p

p = ggplot(subset(dat, channel == "Event_length")[tsne_out_scaled_sub$Y], aes(x=bh_V1, y=bh_V2, color=condition))+
  geom_point(alpha=0.3, size=2)
p



p = ggplot(subset(dat, channel == "Event_length")[tsne_out_scaled_3d$Y], aes(x=bh_V1, y=bh_V2, color=bh_V3))+
  geom_point(alpha=0.3, size=2)
p
tsne_out_scaled_3d
good_channels
marker = good_channels[1]



marker = 'beadDist'
p = ggplot(subset(dat, channel == marker)[tsne_out_scaled_3d$Y], aes(x=bh_V1, y=bh_V2, color=bb.censor_dat(counts)))+
  geom_point(alpha=1, size=1)+
  ggtitle(marker)
p

ggiris <- qplot(Petal.Width, Sepal.Length, data = iris, color = Species)
ggplotly(ggiris)

ggplotly(p)

pdat = subset(dat[tsne_out_scaled_3d$Y], channel == "CXCR4")

# plot conditions
col = rainbow(5)[as.factor(pdat$condition)]
scatterplot3js(pdat$bh_V1,pdat$bh_V2, pdat$bh_V3, size=0.2, color = col , labels =pdat$condition,stroke=NA )
col = bb.map2colormap(bb.censor_dat(pdat$counts, 0.99), cmap=colorRamp(rev(brewer.pal(11,'Spectral'))))
scatterplot3js(pdat$bh_V1,pdat$bh_V2, pdat$bh_V3, size=0.2, color = col , labels =paste(pdat$condition, pdat$counts),stroke=NA)


cond = c('tsne_out_unscaled', 'tsne_out_scaled', 'tsne_out_scaled_sub','tsne_out_unscaled_sub')

for (co in cond){
  
  page = pdf(file.path(fcs_folder,paste('bhSNE',co,'2.pdf',sep='_')),width = 10, height = 7)
  
  #print(paste(tsne_plot$channels, collapse='\n'))
  
  tsne = get(co)
  p = ggplot(subset(dat, !duplicated(id))[tsne$Y], aes(x=bh_V1, y=bh_V2))+
    geom_point(alpha=1, size=0.7, aes(color=condition))+
    scale_colour_brewer(palette="Paired")+
    ggtitle('Conditions')+
    stat_density2d(group=1, color='black')+
    guides(color=guide_legend(override.aes=list(size=5)))
  print(p)
  
  for (chan in unique(dat$channel)){
    plot_dat = subset(dat, channel == chan)[tsne$Y]
    ibreaks = c( quantile(plot_dat$counts,0),
               quantile(plot_dat$counts,0.25),
               quantile(plot_dat$counts,0.5),
               quantile(plot_dat$counts,0.75),
               max(censor_dat(plot_dat$counts)))
    breaks = (ibreaks-min(ibreaks))/max(ibreaks)
    df_b = data.frame(breaks)
    ibreaks = sapply(ibreaks,as.integer)
    df_b$lab= factor(x=1:length(breaks),level=1:length(breaks),
                     labels = c(paste(as.character(ibreaks[1]),'(= min)'),
                                paste(as.character(ibreaks[2]),'(= 25% quantile)'),
                                paste(as.character(ibreaks[3]),'(= median'),
                                paste(as.character(ibreaks[4]),'(= 75% quantile)'),
                                paste(as.character(ibreaks[5]),'(= max)')))
    df_b$x= 0
    df_b$y = 0
    cols =c('darkorchid4','cornflowerblue','grey','yellow','red')
    p = ggplot(plot_dat, aes(x=bh_1, y=bh_2))+
      geom_point(alpha=1, size=0.8,aes( color=censor_dat(counts)))+
      #stat_density2d(group=1, color='black')+
      scale_colour_gradientn(colours=cols,
                             values =breaks)+
      geom_point(data=df_b, aes_string(x='x',y='y',fill='lab'),
                 color='grey',alpha=0, shape=21)+
      scale_fill_manual(values=cols)+
      guides(fill=guide_legend(title = 'Color meaning',override.aes=list(alpha=1,size=5)))+

      ggtitle(chan)
    print(p)
  }
  
  
  
  dev.off()
}


