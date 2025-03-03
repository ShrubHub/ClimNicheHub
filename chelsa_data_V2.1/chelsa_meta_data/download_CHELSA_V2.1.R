library(stringr)

options(timeout=1600)

wget_dist <- fread(file.path("chelsa_data_V2.1","chelsa_meta_data","envidatS3paths.txt"),header = F)

path_url <- wget_dist$V1

path_local <- str_extract(path_url,"CHELSA.*")

path_local <- file.path("chelsa_data_V2.1",path_local)

path_local <- str_replace(path_local,"-","_")

command <- paste0("download.file(path_url[",1:32,"],path_local[",1:32,"],mode='wb')")

for (i in 1:length(command))eval(parse(text = command[i]))


