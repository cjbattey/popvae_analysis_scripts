#alternate hgdp plot
library(ggplot2);library(raster);library(magrittr);library(plyr);library(data.table)
setwd("~/popvae_dev/")
cbPalette <- c("#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette


ld <- fread("out/hgdp/hgdp_latent_coords.txt")
meta <- fread("data/hgdp/hgdp_sample_data.txt")
ld <- merge(ld,meta,by="sampleID")
ld$region <- factor(ld$region,labels=c("Africa","America","Central/South Asia","East Asia","Europe","Middle East","Oceania"),
                    levels=c("AFRICA","AMERICA","CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST","OCEANIA"))
ld <- subset(ld,region!="America")

#normalize latent and geographic coordinates
ld$LD1norm <- (ld$LD1-mean(ld$LD1))/sd(ld$LD1)
ld$LD2norm <- (ld$LD2-mean(ld$LD2))/sd(ld$LD2)
meanlat <- mean(ld$latitude)
sdlat <- sd(ld$latitude)
meanlong <- mean(ld$longitude)
sdlong <- sd(ld$longitude)
ld$longnorm <- (ld$longitude-meanlong)/sdlong
ld$latnorm <- (ld$latitude-meanlat)/sdlat

ld2 <- data.frame(x=c(-ld$LD2norm,ld$longnorm),y=c(ld$LD1norm,ld$latnorm),
                  space=c(rep("latent",nrow(ld)),rep("geographic",nrow(ld))),
                  region=c(as.character(ld$region),as.character(ld$region)))
ld2$region <- factor(ld2$region,levels=c("Africa","America","Central/South Asia","East Asia","Europe","Middle East","Oceania"))

map <- map_data("world")
map$longnorm <- (map$long-meanlong)/sdlong
map$latnorm <- (map$lat-meanlat)/sdlat


pdf("fig/hgdp/hgdp_latent_v_geography_wmap_noAmericas.pdf",useDingbats = F,width=3.5,height=3.75)
p <- ggplot()+
  coord_cartesian(xlim=c(-2,2.2),ylim=c(-5,2.5))+
  theme(axis.title = element_text(size=8),
        legend.position = "bottom",
        legend.background = element_blank())+
  scale_fill_manual(values = cbPalette,name="Region")+
  scale_color_manual(values = cbPalette,guide=F)+
  #scale_color_distiller(palette = "Greens",name="Longitude",guide=F)+
  #scale_fill_distiller(palette = "Greens",name="Longitude")+
  scale_shape_manual(values=c(22,21),name="Space")+
  #xlim(-2,2.2)+
  xlab("Z(x)")+
  ylab("Z(y)")+
  geom_polygon(data=map,aes(x=longnorm,y=latnorm,group=group),lwd=0.15,color="white",fill="grey")+
  geom_segment(data=ld,aes(x=longnorm,y=latnorm,xend=-LD2norm,yend=LD1norm,color=region),lwd=0.25,alpha=0.7)+
  geom_point(data=ld2,aes(x=x,y=y,shape=space,fill=region),color="black",stroke=0.15)+
  guides(fill=guide_legend(override.aes = list(shape=21,size=3),title.position = "top",
                           keyheight = unit(1,"mm"),
                           keywidth = unit(1,"mm"),
                           ncol=2),
         shape=guide_legend(title.position = "top",override.aes=list(size=4,stroke=0.5,
                                                                     fill=c("white","black"),
                                                                     color=c("black","white")),
                            direction = "vertical",keyheight = unit(2,"mm")))
print(p)
dev.off()


