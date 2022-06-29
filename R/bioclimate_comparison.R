library(spThin)
library(raster)
library(tidyverse)
library(gridExtra)

#import occurence points for each species
parent1.data <- read.csv(file = 'occurrence/Lepisorus_annuifrons.csv')
parent2.data <- read.csv(file = 'occurrence/Lepisorus_uchiyamae.csv')
allopolyploid.data <- read.csv(file = 'occurrence/Lepisorus_yamaokae.csv')

# reduce data point for each to avoid spacial autocorrelation using the package 'spthin'
occ.data<-list()
thin.data<-list()

occ.data[[1]]=na.omit(parent1.data)
occ.data[[2]]=na.omit(parent2.data)
occ.data[[3]]=na.omit(allopolyploid.data)

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

##extract climate variables for occurrence points and background points
#read bioclimate data downloaded from worldclim1.4
climdata2.5 <- getData(name = "worldclim",var = "bio", res = 2.5, path = "data/")
env<- list()
for(i in 1:length(thin.data)){
  env[[i]] = na.omit(raster::extract(climdata2.5, thin.data[[i]][[1]]))
}

data=bind_rows(data.frame(env[[1]]),data.frame(env[[2]]), data.frame(env[[3]]), .id = "table_id")

p1=ggplot(data, aes(x=bio1, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio1")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p2=ggplot(data, aes(x=bio2, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio2")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p3=ggplot(data, aes(x=bio3, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio3")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p4=ggplot(data, aes(x=bio4, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio4")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p5=ggplot(data, aes(x=bio5, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio5")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p6=ggplot(data, aes(x=bio6, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio6")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p7=ggplot(data, aes(x=bio7, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio7")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p8=ggplot(data, aes(x=bio8, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio8")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p9=ggplot(data, aes(x=bio9, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio9")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p10=ggplot(data, aes(x=bio10, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio10")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p11=ggplot(data, aes(x=bio11, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio11")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p12=ggplot(data, aes(x=bio12, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio12")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p13=ggplot(data, aes(x=bio13, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio13")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p14=ggplot(data, aes(x=bio14, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio14")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p15=ggplot(data, aes(x=bio15, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio15")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p16=ggplot(data, aes(x=bio16, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio16")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p17=ggplot(data, aes(x=bio17, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio17")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p18=ggplot(data, aes(x=bio18, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio18")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
p19=ggplot(data, aes(x=bio19, color=table_id))+geom_density(aes(color = table_id, linetype = table_id), show.legend = F)+ggtitle("bio19")+theme_minimal()+theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.title.y = element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, nrow=4)


for(i in 1:(length(data[1,])-1)){
  print(i)
  print(pairwise.wilcox.test(data[,i+1], data$table_id, p.adjust.method = 'bonferroni'))
}
  

