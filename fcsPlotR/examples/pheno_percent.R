
## Load the library
source('/home/vitoz/Git/vz-bod/For Stephan/fcsPlotR/fcsPlotR_library.R')

## Set settings

fileDir = '/home/vitoz/imls-bod/for Stephane/'
plotDir = '/home/vitoz/imls-bod/for Stephane/output/'
cluster_col="PhenoGraph K30 AEER24057"
coundition_name_position = c(1,2,3) # used to construct the condition name from the filename
contidion_name_sep ='_'

file_prefix = paste(cluster_col,'1', sep='_') # for the naming of the output files


#### Script
condDict = bb.getInfoFromFileList(list.files(fileDir),contidion_name_sep, coundition_name_position)
dat <- bb.loadConvertMultiFCS(list.files(fileDir),fileDir,condDict)

# prints all the fcs column names
print(colnames(dat))

setnames(dat, cluster_col, 'cluster')
# calculate the counts of each condition per cluster
counts = dat[, list(count=.N), by=list(condition,cluster)]
setkeyv(counts,c('condition', 'cluster'))

counts[,frac_of_condition := count/sum(count),by=condition]
counts[,frac_of_cluster := count/sum(count),by=cluster]

write.csv(counts,file.path(plotDir,paste(file_prefix, 'counts.csv',sep='_')))

counts.
write.csv(counts,file.path(plotDir,paste(file_prefix, 'counts.csv',sep='_')))

# some standard plots

p = ggplot(data = counts,aes(x=factor(cluster), y=count, fill=condition))+
  geom_bar(position='dodge',stat="identity")+
  xlab('Cluster')+
  ggtitle('Counts per Cluster')

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_overview.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=frac_of_cluster, fill=condition))+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Fraction of cluster')+
  ggtitle('Fraction of cluster per cluster')

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_cluster.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=frac_of_condition, fill=condition))+
  geom_bar(position='dodge',stat="identity")+
  xlab('Cluster')+
  ylab('Fraction of condition')+
  ggtitle('Fraction of cluster per condition')

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_condition.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=frac_of_cluster, fill=condition))+
  facet_wrap(~cluster,scale='free_x')+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Fraction of cluster')+
  ggtitle('Fraction of cluster per cluster')

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_cluster_facet.png',sep='_')))

p = ggplot(data = counts,aes(x=condition, y=frac_of_cluster, fill=condition))+
  facet_wrap(~cluster,scale='free_y')+
  geom_bar(position='stack',stat="identity")+
  xlab('Condition')+
  ylab('Fraction of cluster')+
  ggtitle('Fraction of cluster per cluster')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_cluster_facet_dodge.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=frac_of_condition, fill=condition))+
  facet_grid(condition~.,scale='free_y')+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Fraction of Condition')+
  ggtitle('Fraction of Condition per cluster')
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_condition_facet.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=frac_of_condition, fill=condition))+
  facet_grid(condition~cluster,scale='free_x')+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Fraction of Condition')+
  ggtitle('Fraction of Condition per cluster')
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_condition_facet2.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=count))+
  facet_grid(condition~cluster,scale='free')+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Counts')+
  ggtitle('Counts per cluster')
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_counts_of_condition_facet2.png',sep='_')))

p = ggplot(data = counts,aes(x=factor(cluster), y=count))+
  facet_grid(condition~cluster,scale='free_x')+
  geom_bar(position='stack',stat="identity")+
  xlab('Cluster')+
  ylab('Counts')+
  ggtitle('Counts per cluster')
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(plotDir,paste(file_prefix, 'cluster_counts_of_condition_facet3.png',sep='_')))


# p = ggplot(data = counts,aes(x=1, y=frac_of_condition, fill=factor(cluster)))+
#   facet_grid(condition~.,scale='free')+
#   geom_bar(position='stack',stat="identity")+
#   xlab('Cluster')+
#   ylab('Fraction of Condition')+
#   ggtitle('Fraction of Condition per cluster')
# #theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggsave(file.path(plotDir,paste(file_prefix, 'cluster_fractions_of_condition_facet3.png',sep='_')))

p = ggplot(data = counts,aes(x=condition, y=count, fill=condition))+
  facet_wrap(~cluster,scale='free_y')+
  geom_bar(position='dodge',stat="identity")+
  xlab('Cluster')+
  ylab('Counts')+
  ggtitle('Counts per cluster')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(plotDir,paste(file_prefix, 'cluster_counts_facet.png',sep='_')),height=7,width=8)



