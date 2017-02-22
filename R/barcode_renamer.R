if(interactive()){
###
# A script to mass rename fcs files

file_type = '.fcs'

get_substr <- function(x, pos=1, sep='_'){
  return (strsplit(x, sep)[[1]][pos])
}

print('Choose a file in the directories to be renamed.')
fn_directory = dirname(file.choose())

file_list = list.files(path=fn_directory, pattern="\\.fcs$")

print('files found:')
print(file_list)
pos = as.integer(readline(prompt='at which position (sep by _) is the well information?'))
pos_pl = as.integer(readline(prompt='at which position (sep by _) is the plate information?'))


wells = sapply(file_list, function(x){
  x = paste(get_substr(x, pos=c(pos, pos_pl)))
  x = gsub('.fcs','', x)
  return(x)})

print('Veryify well information:')
print(wells)

if (any(duplicated(wells))){
  print('Duplicated wells:')
  print(wells[duplicated(wells)])
  stop('There are duplicated wells in the directory!')


}

print('Choose the csv containing the filename information')
fn_csv = file.choose()

rename_csv = read.csv(fn_csv, header=T)

print('Select the column with the well information:')
well_col = select.list(colnames(rename_csv), preselect = colnames(rename_csv)[1])

print('Select the column with the plate information:')
plate_col = select.list(colnames(rename_csv), preselect = colnames(rename_csv)[2])


print('Select the column with the renamed information:')
rename_col = select.list(colnames(rename_csv), preselect = colnames(rename_csv)[3])

# veryfiy that all wells are in  the csv
rename_dict = as.vector(rename_csv[,rename_col])
names(rename_dict) <- paste(as.vector(rename_csv[,plate_col]),as.vector(rename_csv[,well_col]))
if (!all(wells %in% names(rename_dict) )){
  print('Following wells are missing in the CSV:')
  print(wells[!wells %in% names(rename_dict)])
  stop('Not all wells from the filenames are in the csv!')
}

if (any(duplicated(names(rename_dict) ))){
  stop('There are duplicated wells in the csv files!')
}


fnparts = sapply(strsplit(file_list[1], '_'), function(x) gsub('.fcs','', x))

print('Which positions do you want to keep for the final filename?')
names(fnparts) <- 1:length(fnparts)
pos = as.integer(names(select.list(fnparts, multiple = T)))

print('Do you want to rename the files or save a copy with the new name?')
mode = select.list(c('rename', 'copy'))


do_loop=T
do_rename=F
while (do_loop){
  # rename the files:
  print('old  -  new')
  for (i in 1:length(wells)){
    old_fn = file_list[i]
    new_fn = gsub('.fcs','', old_fn)
    new_fn =paste(strsplit( new_fn, '_')[[1]][pos], collapse='_')
    new_fn = paste(new_fn, rename_dict[wells[i]], sep='_')
    new_fn= paste(new_fn, '.fcs', sep='')

    old_fn = file.path(fn_directory, old_fn)
    new_fn = file.path(fn_directory, new_fn)
    if (do_rename){
      if (mode == 'rename'){
        file.rename(old_fn, new_fn)
      } else {
        file.copy(old_fn, new_fn)
      }
      do_loop = F
    } else {
    print(paste(old_fn, new_fn, sep='  -  '))
    }
  }
  if (do_loop){
  print('Do you want to apply the renaming?')
  do_rename = select.list(c('yes', 'no'))=='yes'
  }
  if (do_rename == 'no'){
    stop('Renaming aborted')
  }
}
}