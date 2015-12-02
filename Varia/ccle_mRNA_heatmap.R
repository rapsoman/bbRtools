library(ggplot2)
library(data.table)
library(gplots)
library(reshape2)
library(tools)


######## Goal ####
# This script makes heatmaps for mRNAs specified in the input_file as well
# as cellines specified bellow. It uses the mRNA data from the Cancer Cell Line Encyclopedia
# (http://www.broadinstitute.org/ccle/data/browseSamples?actionMethod=pages%2Fhome.xhtml%3AbrowseSamplesBean.checkSkipFirstStep%28%29&conversationPropagation=begin)
# as discribed in: http://www.nature.com/nature/journal/v483/n7391/full/nature11003.html
#
# To use the script, SAVE A COPY OF THE SCRIPT and modify things in the COPY
# Usually changing the parameters in the SETUP part should be sufficient.
# Otherwise chaning the heatmap plotting could also be a good idea
#
# The microarray used for had multiple probes per mRNA.
#
# The output will be 4 heatmaps with row normalized and unnormalized values as well as clustered or sorted cellines
#################
##### setup ####

# indicate the output folder
out_dir = '/mnt/imls-bod/LabCode/bbRtools/ExampleData/ExampleOutput'

# Indicate the location of the input file
# this input file needs to be a comma seperated text file with specific columns
input_file = '/mnt/imls-bod/LabCode/bbRtools/ExampleData/150906_mRNAexpression_proteins.csv'


# Look up the correct cell line name here: 
# http://www.broadinstitute.org/ccle/data/browseSamples?actionMethod=pages%2Fhome.xhtml%3AbrowseSamplesBean.checkSkipFirstStep%28%29&conversationPropagation=begin
target_celllines = c('CAKI1', 'ACHN', 'A431', 'Du145', 'PC3')




######## Start script ####
## Load General resources ####
# ccle expression dataset
ccle_filename ='/mnt/imls-bod/Data Vito/Ressources/mRNAdata/CCLE_Expression_2012-09-29.res'
ccle_meta_filename = '/mnt/imls-bod/Data Vito/Ressources/mRNAdata/CCLE_Expression.Arrays.sif_2012-10-18.txt'



####### Start script ########
# determine basename for output
outname =  basename(file_path_sans_ext(input_file))

# load the target protein csv
target_proteins = fread(input_file, stringsAsFactors = F)
target_proteins = melt(target_proteins,id = c('protein','uniprot'))
target_proteins = subset(target_proteins, value != '')

target_proteins_dict = target_proteins[, paste(protein, value, sep= ' - ')]
names(target_proteins_dict) = target_proteins[ , value]
# andrea
#target_celllines = c('LNCaP clone FGC', 'VCaP', 'RWPE1', 'DU 145', 'PC-3')
# stef


# script
cNames <- read.table(ccle_filename, nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE)
#cNames = colnames(fread(ccle_filename, header = T,skip=0,sep='\t',stringsAsFactors = F))
data.expression = fread(ccle_filename, header = F,sep='\t',stringsAsFactors = F)
data.meta = fread(ccle_meta_filename)

# Acctually there is the problem, as there are always two columns (expression level, present or absent)
# per cell line, but only one entry per header
# Replace each second file name in Cnames by name_present and make a logical vector out of the peresent P absent a

tmp= paste('present',cNames[2+seq(1,length(cNames)-2,2)],sep='_')
cNames[3+seq(1,length(cNames)-2,2)] = tmp


setnames(data.expression,colnames(data.expression),as.matrix(cNames))

### Match the lab cellines with the database ####
data.meta$'Cell line primary name'[data.meta$'Cell line primary name' %in% target_celllines]
data.meta[, matchName := tolower(get('Cell line primary name'))]
data.meta[, matchName := gsub('-|\ ','',matchName)]

target_celllines.matchName = gsub('-|\ ','',tolower(target_celllines))

data.meta[, inLab := F]
data.meta[ matchName %in% target_celllines.matchName, inLab := TRUE]

#subset the data expression to include only the target cell lines
targetCCLE = subset(data.meta,inLab == TRUE)
targetCCLE = targetCCLE[,get('CCLE name')]

data.expression = data.expression[,colnames(data.expression) %in% c('Description',targetCCLE),with=F]

### Match the protein of interest with transcript of interest ####



setkey(data.expression,Description)
unique(data.expression[ Description %in% target_proteins$value,Description])

#subset the data expression to include only the target proteins
data.expression = subset(data.expression, Description %in% target_proteins$value)



### Make a clustermap to find look at the variability between the cell lines ####

rowNam =data.expression[,Description]
dataMat = data.matrix(data.expression,rownames.force=T)[,-1]
rownames(dataMat)=target_proteins_dict[rowNam]
dataMat = dataMat[, order(colnames(dataMat))]
whclust <- function(x){hclust(x, method='ward')}

pdf(file.path(out_dir, paste(outname,'CellllinesVsMRNA_unnorm_clustered.pdf',sep='_')),
    height=20, width=10)
heatmap.2((dataMat),
          scale ='none',
          trace = "none",
          col=greenred(75),
          hclustfun =whclust,
          density.info ='none',
          Rowv=T,
          #Colv =NA,
          #keyorient=2,
          xlab = 'Marker',
          sepwidth = c(0,0),
          cexCol=0.6,
          cexRow = 0.8,
          offsetRow = 0.1,
          offsetCol = 0.1,
          margins = c(10,20) ,
          lhei=c(.1,0.9)
          #ylab ='Genes',
          #comments = data.frame(names = row.names(tab_Prot))
          
)

dev.off()

pdf(file.path(out_dir, paste(outname,'CellllinesVsMRNA_rownorm_clustered.pdf',sep='_')),
    height=20, width=10)
heatmap.2((dataMat),
          scale ='row',
          trace = "none",
          col=greenred(75),
          hclustfun =whclust,
          density.info ='none',
          Rowv=T,
          #Colv =NA,
          #keyorient=2,
          xlab = 'Marker',
          sepwidth = c(0,0),
          cexCol=0.6,
          cexRow = 0.8,
          offsetRow = 0.1,
          offsetCol = 0.1,
          margins = c(10,20) ,
          lhei=c(.1,0.9)
          #ylab ='Genes',
          #comments = data.frame(names = row.names(tab_Prot))
          
)

dev.off()

pdf(file.path(out_dir, paste(outname,'CellllinesVsMRNA_rownorm_sorted.pdf',sep='_')),
    height=20, width=10)

heatmap.2((dataMat),
          scale ='row',
          trace = "none",
          col=greenred(75),
          hclustfun =whclust,
          density.info ='none',
          Rowv=F,
          Colv =T,
          #keyorient=2,
          xlab = 'Marker',
          sepwidth = c(0,0),
          cexCol=0.6,
          cexRow = 0.8,
          offsetRow = 0.1,
          offsetCol = 0.1,
          margins = c(10,20) ,
          lhei=c(.1,0.9)
          #ylab ='Genes',
          #comments = data.frame(names = row.names(tab_Prot))
          
)

dev.off()

pdf(file.path(out_dir, paste(outname,'CellllinesVsMRNA_unnorm_sorted.pdf',sep='_')),
    height=20, width=10)
heatmap.2((dataMat),
          scale ='none',
          trace = "none",
          col=greenred(75),
          hclustfun =whclust,
          density.info ='none',
          Rowv=F,
          Colv =T,
          #keyorient=2,
          xlab = 'Marker',
          sepwidth = c(0,0),
          cexCol=0.6,
          cexRow = 0.8,
          offsetRow = 0.1,
          offsetCol = 0.1,
          margins = c(10,20) ,
          lhei=c(.1,0.9)
          #ylab ='Genes',
          #comments = data.frame(names = row.names(tab_Prot))
          
)
dev.off()


# Write out the info for the cell lines as a csv
#meltDat = melt(data.expression,id.vars='variable')
#setnames(meltDat,c('variable','variable.1'),c('condition','channel'))

#write.csv(meltDat,file=file.path(out_dir, paste(outname,'CellllinesVsMRNA_rawdata.txt',sep='_')))
