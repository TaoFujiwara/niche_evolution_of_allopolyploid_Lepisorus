library(sp)
library(raster)
library(maptools)
library(rgdal)
library(dismo)
library(spThin)
library(ecospat)
library(gpclib)
library(rgeos)
library(sf)
library(ecospat)
library(tidyverse)


## preparation of occurrence data
occ.data <- list()
thin.data <- list()

# read raw data of occurrence data for each species
occ.data[[1]]=na.omit(read.csv('occurrence/Lepisorus_annuifrons.csv'))
occ.data[[2]]=na.omit(read.csv('occurrence/Lepisorus_uchiyamae.csv'))
occ.data[[3]]=na.omit(read.csv('occurrence/Lepisorus_yamaokae.csv'))

#thinning occurence data to reduce an effect of spatial autocorrelation
for (i in 1:length(occ.data)){
  thin.data[[i]] = spThin::thin(occ.data[[i]], verbose = FALSE, 
                        lat.col = "decimalLatitude",
                        long.col = "decimalLongitude",
                        spec.col = "species",
                        thin.par=10,
                        reps = 1, 
                        write.files = FALSE,
                        write.log.file = FALSE,
                        locs.thinned.list.return = TRUE)
}

##Preparation of background regions
#shp file of world map delimited by countries was downloaded from  Natural Earth (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/)
wrld=sf::read_sf("map/ne_10m_admin_0_countries.shp")
#extract Japan from world map
jpn<-wrld%>%
dplyr::filter(ADM0_A3 %in% "JPN")
#plot to check the polygon
ggplot(data=jpn)+geom_sf()+theme_void()


#read shp vector file of SDM for each species, converted from .asc raster file in QGIS  
sdm.shp <- list()
#1,2,3 are L. annuifrons, L. ushiyamae and L. yamaokae respectively
sdm.shp[[1]]<-sf::read_sf("shp/Lepisorus_annuifrons.shp")
sdm.shp[[2]]<-sf::read_sf("shp/version1/uchiyamae.shp")
sdm.shp[[3]]<-sf::read_sf("shp/version1/yamaokae.shp")

#make new shp file of shared area between Japan contour and each species's SDM. 
back <- list()
for (i in 1:length(sdm.shp)){
  st_crs(sdm.shp[[i]]) = "+proj=longlat +datum=WGS84"
  temp=st_make_valid(sdm.shp[[i]])
  back[[i]]=st_intersection(st_union(temp),jpn$geometry)%>% #clip shared areas between SDM and Japan polygons
            st_cast("POLYGON")%>% #MULTIPOLYGON to POLYGON
            as_Spatial() #covert sf object to spatial object
  plot(back[[i]])
}

#sampling 1000 random points within background region for each species
back.data<-list()
n=1000
for (i in 1:length(back)){
  back.data[[i]]=spsample(back[[i]], n=n, "random")
}

##extract climate variables for occurrence points and background points
#read bioclimate data downloaded from worldclim1.4
climdata2.5 <- getData(name = "worldclim",var = "bio", res = 2.5, path = "data/")
reduced_climdata2.5 <- getData(name = "worldclim",var = "bio", res = 2.5, path = "data/")

env.back <- list()
env.real <- list()
for(i in 1:length(thin.data)){
  env.real[[i]] = na.omit(raster::extract(climdata2.5, thin.data[[i]][[1]]))
  env.back[[i]] = na.omit(raster::extract(climdata2.5, back.data[[i]]))
}

#make output folder for whole outputs
whole_output_dir="test_output"
dir.create(path = whole_output_dir)

all_output_dir="all"
path_to_all = paste(whole_output_dir,all_output_dir, sep="/")
dir.create(path = path_to_all)

#Put all datasets together one list
all.env.back <- do.call(rbind.data.frame, env.back)
all.env.real <- do.call(rbind.data.frame, env.real)
data.env <- rbind(all.env.real, all.env.back)

#make labels to distinguish between real occurrences and background points
w <- c(rep(0, nrow(all.env.real)), rep(1, nrow(all.env.back)))

#PCA
pca.cal <- dudi.pca(data.env, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
pdf_name=paste(path_to_all, "ecospat.plot.contrib.pdf", sep="/")
pdf(pdf_name)
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)
dev.off()

#Display what bioclimate variable contribute to the clustering
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

#Sort each score in PCA into each of occurrence and background points for each species. 
list <- c(rep(0, nrow(env.real[[1]])), rep(1, nrow(env.real[[2]])), rep(2, nrow(env.real[[3]])), rep(3, nrow(env.back[[1]])),rep(4, nrow(env.back[[2]])), rep(5, nrow(env.back[[3]])))
  
sort <- list()
sep.data <- list()
  
for(i in 1:length(as.numeric(levels(as.factor(list))))){
  temp=c()
  for(j in 1:length(list)){
    if(list[j]==i-1){
      temp=append(temp,j)
      sort[[i]]<-temp
    }
  }
}
  
for(i in 1:length(sort)){
  sep.data[[i]]=pca.cal$li[sort[[i]],]
}
  
sp.pca.data=list()
back.pca.data=list()
  
sp.pca.data[[1]]=sep.data[[1]]
sp.pca.data[[2]]=sep.data[[2]]
sp.pca.data[[3]]=sep.data[[3]]
sp.pca.data[[4]]=rbind(sp.pca.data[[1]],sp.pca.data[[2]])
back.pca.data[[1]]=sep.data[[4]]
back.pca.data[[2]]=sep.data[[5]]
back.pca.data[[3]]=sep.data[[6]]
back.pca.data[[4]]=rbind(sep.data[[4]][1:(length(sep.data[[4]][,1])),],sep.data[[5]][1:(length(sep.data[[5]][,1])),])
total.scores.back <- do.call(rbind.data.frame, sep.data)
  
R <- 1000
z <- list()
for(i in 1:length(sp.pca.data)){
  z[[i]] <- ecospat.grid.clim.dyn(glob=total.scores.back, glob1=back.pca.data[[i]], sp=sp.pca.data[[i]], R = R, kernel.method='ks')
}

#Plot niche space for each species
for(i in 1:length(z)){
  ecospat.plot.niche(z[[i]])
  pdf_name=paste(path_to_all,paste('ecospat.plot.niche',as.character(i), ".pdf",sep=""),sep="/")
  pdf(pdf_name)
  ecospat.plot.niche(z[[i]])
  dev.off()
}

##niche overlap
#preparing a list for pairs to compare
#[[1]]:L. annuifrons v.s. L.uchiyamae
#[[2]]:L. annuifrons v.s. L. yamaokae
#[[3]]:L. uchiyamae v.s. L. yamaokae
#[[4]]:expected additive niche v.s. realized niche in allopolyploid

pairs <- list(c(1,2),c(1,3),c(2,3),c(3,4))

for(i in 1:length(pairs)){
  #Niche overlap
  D.overlap=ecospat.niche.overlap(z[[pairs[[i]][1]]], z[[pairs[[i]][2]]], cor=T)
  output=paste(path_to_all, paste("D.overlap",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".csv", sep=""), sep="/")
  write.csv(D.overlap, output)
  
  #Perform the Niche Equivalency Test with ecospat.niche.equivalency.test() according to Warren et al. (2008)
  #Niche equivalency test H1: Is the overlap between two species niche lower than two random niches?
  eq.test=ecospat.niche.equivalency.test(z[[pairs[[i]][1]]], z[[pairs[[i]][2]]], rep=1000, alternative = "lower", ncore=8)
  output=paste(path_to_all, paste("eq.test",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".csv", sep=""), sep="/")
  write.csv(eq.test, output)
  
  ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
  pdf_name=paste(path_to_all, paste("'ecospat.niche.equivalency.test'",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".pdf", sep=""),sep="/")
  pdf(pdf_name)
  ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
  dev.off()
  
  sim.lower.test=ecospat.niche.similarity.test(z[[pairs[[i]][1]]], z[[pairs[[i]][2]]], rep=1000, alternative = "lower", ncore=8, rand.type=1)
  output=paste(path_to_all, paste("sim.lower.test",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".csv", sep=""), sep="/")
  write.csv(sim.lower.test, output)
  
  ecospat.plot.overlap.test(sim.lower.test, "D", "Similarity")
  pdf_name=paste(path_to_all, paste("'sim.lower.test'",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".pdf", sep=""),sep="/")
  pdf(pdf_name)
  ecospat.plot.overlap.test(sim.lower.test, "D", "Similarity")
  dev.off()
  
  sim.greater.test=ecospat.niche.similarity.test(z[[pairs[[i]][1]]], z[[pairs[[i]][2]]], rep=1000, alternative = "greater", ncore=8, rand.type=1)
  output=paste(path_to_all, paste("sim.greater.test",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".csv", sep=""), sep="/")
  write.csv(sim.greater.test, output)
  
  ecospat.plot.overlap.test(sim.greater.test, "D", "Similarity")
  pdf_name=paste(path_to_all, paste("sim.greater.test",as.character(pairs[[i]][1]),as.character(pairs[[i]][2]),".pdf", sep=""),sep="/")
  pdf(pdf_name)
  ecospat.plot.overlap.test(sim.greater.test, "D", "Similarity")
  dev.off()
  
}

#exploration of niche evolution of allopolyploid 
#by comparing between expected and realized niches under USE framework

ecospat.plot.niche.dyn(z[[4]],z[[3]], quant=0.25, interest=2,title= "Niche Overlap", name.axis1="PC1",name.axis2="PC2")
ecospat.niche.dyn.index(z[[3]],z[[4]])
pdf_name=paste(path_to_all, paste("USE",".pdf", sep=""),sep="/")
pdf(pdf_name)
ecospat.plot.niche.dyn(z[[4]],z[[3]], quant=0.25, interest=2,title= "Niche Overlap", name.axis1="PC1",name.axis2="PC2")
dev.off()







