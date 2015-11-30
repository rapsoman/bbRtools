library(igraph) # for graph
library(data.table) # for fast data tables
library(ggplot2) # For plotting
library(rjson) # For reading json
library(RColorBrewer) # For nice color palettes
library(LSD)
library(colorRamps)
library(fields)

net_file = '/mnt/imls-bod/Daniel_Data/Macrophage network/Network_cyto.cyjs'
data_file = '/mnt/imls-bod/Daniel_Data/Macrophage network/20151120 Signaling 1.04.csv'
out_folder = '/mnt/imls-bod/Daniel_Data/Macrophage network/'

#load json
jdat = fromJSON(file = net_file)


load_json_graph = function(jdat){
  "
  Loads a graph with metadata from a cytoscape.js JSON file
  "
  edgelist = sapply(jdat$elements$edges, function(e){
    # Loads the edge matrix f  
    return(c(e$data$source, e$data$target))
  })
  # Create graph
  g = graph_from_edgelist(t(edgelist))
  
  # Loop over the nodes to get metadata
  for (v in jdat$elements$nodes){
    
    idx = which(V(g)$name == v$data$id)
    
    # Loop over the metadata entries
    for (at in names(v$data)){
      
      g = set.vertex.attribute(g, at, idx, v$data[[at]])
    }
    
    # set position metadata
    g = set.vertex.attribute(g,'posx',idx, as.numeric(v$position$x))
    g = set.vertex.attribute(g,'posy',idx, as.numeric(v$position$y))
    
  }
  
  # Loop over edges metadata
  for (idx in 1:length(jdat$elements$edges)){
    e = jdat$elements$edges[[idx]]
    for (at in names(e$data)){
      g = set.edge.attribute(g,at,idx, e$data[[at]])
    }
  }
  return (g)
}

g = load_json_graph(jdat)


list.edge.attributes(g)
list.vertex.attributes(g)

plot.igraph(g,
            layout=-cbind(V(g)$posx, V(g)$posy),
            rescale=T,
            edge.arrow.size=0.1,
            vertex.size =7,
            vertex.label.cex = 0.5)

# Make a function to map interaction value to color
arrow_col_map = function(val){
  val = as.numeric(val)
  if (is.na(val)){
    val = 1
  }
  if (val < 0){
    return('red')
  } else {
    return('black')
  }
}

arrow_lty_map = function(val){
  val = abs(as.numeric(val))
  if (is.na(val)){
    val = 1
  }
  if (val == 3){
    return(3)
  } else {
    return(1)
  }
}


# Color the shit:
display.brewer.all()
# look up http://www.r-bloggers.com/r-using-rcolorbrewer-to-colour-your-figures-in-r/
# brewer.pal(8,"Set3")
plot.igraph(g,
            layout=-cbind(V(g)$posx, V(g)$posy),
            rescale=T,
            edge.arrow.size=0.1,
            vertex.size =7,
            vertex.shape = 'square',
            vertex.label.cex = 0.5,
            edge.color = sapply(E(g)$interaction, arrow_col_map),
            edge.lty=sapply(E(g)$interaction, arrow_lty_map))

E(g)[is.na(E(g)$interaction)]
E(g)$interaction2 =2

E(g)[is.na(E(g)$interaction)]


### Map data ####
# read data to be mapped
map_dat = fread(data_file,stringsAsFactors = F)

# set the row name
setnames(map_dat, 'V1', 'antibody')

# set the row name as key
setkey(map_dat, 'antibody')

# # ref condition names
# ref_conditions = c("MCSF1.fcs - MCSF1.fcs" ,              "MCSF2.fcs - MCSF2.fcs"         ,      "MCSF3.fcs - MCSF3.fcs"    )

# 
map_dat = melt.data.table(map_dat, id.vars = 'antibody', variable.name = 'condition', value.name = 'counts', variable.factor = F)

# get reference condition as a logic
map_dat[, ref_filter := grepl('ctrlMCSF', condition)]

# divide by the mean of all reference conditions
map_dat[ ,  foldchange_toref := counts/ mean(counts[ref_filter]), by=antibody]

map_dat[ , log2_foldchange_toref := log2(foldchange_toref)]

## split conditions and timepoints
map_dat[, stimulation:= sapply(condition, function(x){
  return(strsplit(x, '_')[[1]][1])
})]

map_dat[, timepoint := sapply(condition, function(x){
  return(as.numeric(strsplit(x, '_')[[1]][2]))
})]

# get reference condition as a logic

map_dat[ , norm_per_tp := log2(counts/mean(counts)) , by=list(antibody, timepoint) ]

map_dat[ , z_per_tp := (counts-mean(counts))/sd(counts) , by=list(antibody, timepoint) ]




# find antibodies in the network
map_dat[antibody %in% V(g)$FCS_name, antibody] 

# set double key with a list of characters of the column names
setkey(map_dat, antibody, condition)

# entries can now be accessed by the keys as follows
map_dat[.('p90RSK_phospho_230 - Panel 1','GC-TGFb_15_min'),]

# without indexing this would have looked like
map_dat[ antibody =='Erk1_2_phospho_105 - Panel 1' & condition =='MCSF1.fcs - MCSF1.fcs',  ]
map_dat[ condition =='GC-TGFb_15_min',  ]
map_dat[counts < 2,]
map_dat[ , condition]

#### Map the values of the condition to the network ####

V(g)$value = NA
cur_condition = condition_names[1]

condition_names = map_dat[, unique(condition[order(stimulation, timepoint)])]
# loop over all conditions and print to pdf


# color mapping
vertex_color_map = function(x, xmax=NaN, symmetric=NA, cmap=colorRamp(c('blue','white', 'red')), na_col='cornsilk3'){
  if (is.na(xmax)){
    x_max = max(abs(x), na.rm = T)
  }
  
  x = x/xmax
  
  # if any value is negative, make symmetric by scaling
  if ((!is.na(symmetric) & symmetric == T) | any(x[!is.na(x)] <0)){
    x = (x+1)/2
  }
  
  # apply
  c_x =cmap(x)
  
  # deal with NA
  if (is.character(na_col)){
    na_rgb_col = col2rgb(na_col)
  } else {
    na_rgb_col = na_col
  }
  
  if (any(is.na(c_x))){
    na_fil = which(apply(is.na(c_x),1, any))
    
    for (row in na_fil){
      c_x[row,] = na_rgb_col
    }
  }
  
  rgb_x = rgb(c_x/255)
  return(rgb_x)
}

##### 
# Settings for plotting
display.brewer.all()
value_var = 'log2_foldchange_toref'
cmap = colorRamp(rev(brewer.pal(11,'Spectral')))
#cmap = colorRamp(c('blue','lightblue', 'white', 'tomato','red'))
# calculate the max fold change overall for all antibodies
xmax = map_dat[, max(abs(get(value_var)))]


setkeyv(map_dat, c('antibody', 'stimulation','timepoint'))
####
pdf(file = file.path(out_folder,paste('networks_mapped', value_var, 'log2_foldchange_toref.pdf', sep='_')), width=20, height=10)
for (stim in map_dat[, unique(stimulation)]){
  
  # get the unique timepoints
  unitp = map_dat[stimulation %in% stim, sort(unique(timepoint))]
  
  # initialize figure pannel
  par(mfrow=c(1,length(unitp)))
  
  
  # loop through and plot all timepoints
  for (tp in unitp){
    
    
    
    #### Do the plotting #### 
    cur_condition = paste(stim, tp, sep='_')
    
    
    # Add the current measurement to all nodes. A "." may be used to call keys instead of "list"
    for (idx in 1:length(V(g))){
      cur_v = V(g)[[idx]]
      #if the Id in the graph data matches a measurement then we fill the graph data "g" with "value_var"
      if (cur_v$FCS_name %in% map_dat$antibody){
        g = set.vertex.attribute(g,'value',idx, map_dat[.(cur_v$FCS_name, stim , tp),get(value_var)])
      }
    }
    

    
    p = plot.igraph(g,
                    layout=-cbind(V(g)$posx, V(g)$posy),
                    rescale=T,
                    edge.arrow.size=0.2,
                    vertex.size =7,
                    vertex.shape = 'square',
                    vertex.label.cex = 0.5,
                    edge.color = sapply(E(g)$interaction, arrow_col_map),
                    edge.lty=sapply(E(g)$interaction, arrow_lty_map),
                    vertex.color =vertex_color_map(V(g)$value, symmetric=T, xmax=xmax, cmap=cmap))
    print(p)
    title(cur_condition)
    image.plot( legend.only=TRUE, zlim=c(-xmax,xmax), col=vertex_color_map(1:200, xmax=200, cmap=cmap), legend.shrink = 0.2) 
    
  }}
dev.off()









# #####
# pdf(file = file.path(out_folder,paste('networks_mapped', value_var, '.pdf', sep='_')), width=15, height=15)
# for (cur_condition in condition_names){
#   for (idx in 1:length(V(g))){
#     cur_v = V(g)[[idx]]
#     if (cur_v$FCS_name %in% map_dat$antibody){
#       g = set.vertex.attribute(g,'value',idx, map_dat[list(cur_v$FCS_name, cur_condition),get(value_var)])
#     }
#   }
#   
#   # color mapping
#   
#   
#   
#   vertex_color_map = function(x, xmax=NaN, symmetric=NA, cmap=colorRamp(c('blue','white', 'red')), na_col='cornsilk3'){
#     if (is.na(xmax)){
#       x_max = max(abs(x), na.rm = T)
#     }
#     
#     x = x/xmax
#     
#     # if any value is negative, make symetric by scaling
#     if ((!is.na(symmetric) & symmetric == T) | any(x[!is.na(x)] <0)){
#       x = (x+1)/2
#     }
#     
#     # apply
#     c_x =cmap(x)
#     
#     # deal with NA
#     if (is.character(na_col)){
#       na_rgb_col = col2rgb(na_col)
#     } else {
#       na_rgb_col = na_col
#     }
#     
#     if (any(is.na(c_x))){
#       na_fil = which(apply(is.na(c_x),1, any))
#       
#       for (row in na_fil){
#         c_x[row,] = na_rgb_col
#       }
#     }
#     
#     rgb_x = rgb(c_x/255)
#     return(rgb_x)
#   }
#   
#   
#   p = plot.igraph(g,
#                   layout=-cbind(V(g)$posx, V(g)$posy),
#                   rescale=T,
#                   edge.arrow.size=0.2,
#                   vertex.size =7,
#                   vertex.shape = 'square',
#                   vertex.label.cex = 0.5,
#                   edge.color = sapply(E(g)$interaction, arrow_col_map),
#                   edge.lty=sapply(E(g)$interaction, arrow_lty_map),
#                   vertex.color =vertex_color_map(V(g)$value, symmetric=T, xmax=xmax, cmap=cmap))
#   title(cur_condition)
#   image.plot( legend.only=TRUE, zlim=c(-xmax,xmax), col=vertex_color_map(1:200, xmax=200, cmap=cmap), legend.shrink = 0.5) 
#   print(p)
# }
# dev.off()