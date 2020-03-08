# install.packages("rJava")
# library(devtools)
# install_github("johnbaums/rmaxent")
# library("rmaxent")
# get_maxent(version = "latest", quiet = FALSE)
# install_github("danlwarren/ENMTools")
## If installation doesn't go, install ecospat

## Set working directory, which must contain the folder with asc files
wdpath <- "/Users/oscar/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_MaxEntModels_Feb042020"
# wdpath <- "C:/Users/OSCAREOSPINATOBON/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_MaxEntModels_Feb042020"
setwd(wdpath)

## Set directory for results without final slash
outpath <- "/Users/oscar/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_NicheOverlap_Feb102020"
# outpath <- "C:/Users/OSCAREOSPINATOBON/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_NicheOverlap_Feb102020"

## Specifiy layers to be excluded (e.g categorical variables)
# exclude <- c("landfire", "soil")

## Specify RegEx for variable names (e.g."bio[0-9]+", "elev"...)
# var.tokens <- c("bio_aligned_[0-9]+", "elev")
var.tokens <- c("bio_aligned_[0-9]+", "elev_aligned_20", "landfireEVC_aligned_21", "landfireEVT_aligned_22", "soil_aligned_23")

## Specify RegEx for two coordinate files
coordfiles.tokens <- c("records_Acaeca", "records_Axera")

## Specify names of two taxa to be tested, which should be part of coordfiles.tokens
species.name <- c("caeca", "xera")

mx.args <- c("togglelayertype=landfireEVC_aligned_21", "togglelayertype=landfireEVT_aligned_22", "togglelayertype=soil_aligned_23", "betamultiplier=2", 
"randomtestpoints=25", "randomseed=TRUE", "replicatetype=bootstrap") ## Parameters for MaxEnt. Note the same parameters will be applied to both species' models

library("ENMTools")
library("stringr")

env.files <- list.files(path = paste(wdpath, "/DATASET_ascii_noncorr", sep=""), pattern = "asc$", full.names = TRUE)
# env.files <- env.files[!grepl(paste(exclude, collapse='|'), env.files)]
env <- stack(env.files)

var.names <- str_match(env.files, paste(var.tokens, collapse='|'))[,1]
names(env) <- var.names
env <- setMinMax(env)

## Look for coordinate files using species.names
coord.files <- list.files(path = wdpath, pattern = "csv$", full.names = TRUE)
coord.files <- coord.files[grepl(paste(coordfiles.tokens, collapse='|'), coord.files)]


## Create species objects for each taxon
species_1 <- enmtools.species(species.name = species.name[1])
species_1$presence.points <- read.csv(grep(species.name[1], coord.files, value=T))[,2:3]
species_1$range <- background.raster.buffer(species_1$presence.points, 50000, mask = env)
species_1$background.points <- background.points.buffer(points = species_1$presence.points, radius = 20000, n = 1000, mask = env[[1]])
crs(species_1$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

species_2 <- enmtools.species(species.name = species.name[2])
species_2$presence.points <- read.csv(grep(species.name[2], coord.files, value=T))[,2:3]
species_2$range <- background.raster.buffer(species_2$presence.points, 50000, mask = env)
species_2$background.points <- background.points.buffer(points = species_2$presence.points, radius = 20000, n = 1000, mask = env[[1]])
crs(species_2$range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

bg.mx.asym_sp1 <- background.test(species.1 = species_1, species.2 = species_2, env = env, type = "mx", nreps = 100, test.type = "asymmetric", low.memory=T, args=mx.args, 
	rep.dir="/Users/oscar/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_NicheOverlap_Feb102020/caeca_vs_xera_Feb102020") 

bg.mx.asym_sp2 <- background.test(species.1 = species_2, species.2 = species_1, env = env, type = "mx", nreps = 100, test.type = "asymmetric", low.memory=T, args=mx.args,
	rep.dir="/Users/oscar/Dropbox/TESIS_AMPHISBAENA/HabitatModels_Amphisbaena_MaxEnt_Feb2020/amphisbaena_NicheOverlap_Feb102020/xera_vs_caeca_Feb102020")

write.table(bg.mx.asym_sp1$reps.overlap, file=paste(outpath,"/",species.name[1],"VS",species.name[2],"_",species.name[1],"_background_overlap_100reps.txt", sep=""), quote=F, sep=",")
write.table(bg.mx.asym_sp1$p.values, file=paste(outpath,"/",species.name[1],"VS",species.name[2],"_",species.name[1],"background_overlap_100reps_PVal.txt", sep=""), quote=F, sep=",")
write.table(bg.mx.asym_sp2$reps.overlap, file=paste(outpath,"/",species.name[1],"VS",species.name[2],"_",species.name[2],"_background_overlap_100reps.txt", sep=""), quote=F, sep=",")
write.table(bg.mx.asym_sp2$p.values, file=paste(outpath,"/",species.name[1],"VS",species.name[2],"_",species.name[2],"background_overlap_100reps_PVal.txt", sep=""), quote=F, sep=",")
