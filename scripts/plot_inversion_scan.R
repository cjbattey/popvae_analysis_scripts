library(ggplot2);library(plyr);library(data.table)
library(plyr);library(sp);library(magrittr);library(MASS)
setwd("~/popvae/")

#load data
sd <- fread("data/ag1000g/anopheles_samples_sp.txt")
files <- list.files("out/ag1000g/inversion_windows",full.names = T)
files <- files[grepl("latent_coords.txt",files)]
for(i in 1:length(files)){
  if(i==1){
    vae <- fread(files[i],header=T)
    vae$start <- basename(files[i]) %>% strsplit("_") %>% unlist() %>% .[1] %>% strsplit("-") %>% unlist() %>% .[1]
    names(vae)[1:2] <- c("LD1","LD2")
  } else {
    tmp <- fread(files[i],header=T)
    tmp$start <- basename(files[i]) %>% strsplit("_") %>% unlist() %>% .[1] %>% strsplit("-") %>% unlist() %>% .[1]
    names(tmp)[1:2] <- c("LD1","LD2")
    vae <- rbind(vae,tmp)
  }
}
vae <- merge(vae,sd,by="sampleID")
vae$species[is.na(vae$species)] <- "unknown"
vae$window <- as.numeric(vae$start)
unique(vae$country)
vae$country[vae$country=="France"] <- "Mayotte"
vae$country[vae$country=="Gambia, The"] <- "The Gambia"
vae$species[vae$species=="coluzzi"] <- "coluzzii"

#scale latent coordinates and get distance matrices
dists <- list()
j <- 1
windows <- unique(vae$start)
for(i in windows){
  a <- as.matrix(subset(vae,start==i)[,c("LD1","LD2")])
  a <- apply(a,2,function(e){
    e <- e-min(e)
    e <- e/max(e)
    return(e)
  })
  b <- as.matrix(dist(a)) #or should we do the covariance matrix from Li & Ralph 2019? 
  dists[[j]] <- b
  j <- j+1
}
names(dists) <- windows

#get distance among distance matrices
nwindows <- length(windows)
pairs <- combn(1:nwindows,2)
dmat <- matrix(nrow=nwindows,ncol=nwindows)
for(i in 1:ncol(pairs)){
  e <- pairs[,i]
  m1 <- dists[[ e[1] ]]
  m2 <- dists[[ e[2] ]]
  d <- sum((m1-m2)^2)
  dmat[e[1],e[2]] <- d
}
dmat <- as.dist(t(dmat))

#run MDS
mds <- cmdscale(dmat,k=1,add=T)
pd <- data.frame(mds=mds$points,window=as.numeric(unique(vae$start)),stringsAsFactors = F)

##################### plots ######################
theme_set(theme_classic()+theme(axis.title=element_text(size=7),
                                axis.text=element_text(size=6),
                                legend.text=element_text(size=6),
                                legend.title=element_text(size=7,face="bold"),
                                strip.text = element_text(size=8),
                                title=element_text(size=8),
                                axis.line = element_line(size=0.35),
                                strip.background = element_blank()))
pal <- c(RColorBrewer::brewer.pal(12,"Paired"),"black","grey","violet","navyblue")
windows_to_plot <- c(2e7,2.06e7,4.2e7,4.24e7)

#color by country
p1a <- ggplot(data=subset(vae,window %in% windows_to_plot),
             aes(x=LD1,y=LD2,fill=country,shape=species))+
  theme(legend.box = "horizontal",
        legend.margin = margin(0,30,-5,0),
        legend.position = "top",
        strip.text = element_blank(),
        plot.background = element_blank())+
  facet_wrap(~start,scales="free",nrow=1)+
  scale_fill_manual(values=pal,name="Country")+
  scale_shape_manual(values=c(21,22,23),name="Species",labels=c(expression(italic(A.~coluzzi)),
                                                                expression(italic(A.~gambiae)),
                                                                "unknown"))+
  geom_point(stroke=0.01,size=1.35)+
  guides(fill=guide_legend(override.aes = list(size=2,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "left",
                           ncol=4),
         shape=guide_legend(override.aes = list(size=2,stroke=0.5),
                            keyheight = unit(1,"mm"),
                            title.position = "left",
                            ncol=1))

#color by species
p1b <- ggplot(data=subset(vae,window %in% windows_to_plot),
             aes(x=LD1,y=LD2,fill=species))+
  theme(legend.box = "horizontal",
        legend.margin = margin(0,0,-7,0),
        legend.position = "top",
        axis.title.x=element_text(vjust=3),
        strip.text = element_blank(),
        plot.background = element_blank())+
  facet_wrap(~start,scales="free",nrow=1)+
  scale_fill_manual(values=pal,name="Species",labels=c(expression(italic(A.~coluzzii)),
                                                       expression(italic(A.~gambiae)),
                                                       "unknown"))+
  geom_point(shape=21,stroke=0.01,size=1.45)+
  guides(fill=guide_legend(override.aes = list(size=3,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "left",
                           ncol=2))


arrows <- data.frame(x=windows_to_plot+1e5)
arrows$y <- pd$mds[pd$window %in% windows_to_plot]
arrows$xend <- c(1.4e7,2.2e7,3.8e7,4.7e7)
arrows$yend <- rep(1e5,4)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

p2 <- ggplot(data=pd,aes(x=window,y=mds))+
  annotate(geom="rect",xmin=20528089,xmax=42165182,ymin=-45000,ymax=-38000,
           color=NA,fill="grey",alpha=0.8)+
  annotate(geom="text",x=2.4e7,y=-30000,label="2La inversion",size=2.6)+
  ylab("MDS Coordinate")+xlab("Chromsome 2L Position (bp)")+
  geom_segment(data=arrows,aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(2,"mm"),type="closed",angle=25),
               lwd=0.35,color="grey30")+
  geom_step(lwd=0.4)+
  scale_x_continuous(label=scientific,minor_breaks = seq(1e7,5e7,1e7))

library(cowplot)
pdf("~/popvae/fig/ag1000g_inversion_scan_by_country.pdf",width=6.5,height=3.5,useDingbats = F)
ggdraw()+
  draw_plot(p2,0,0,1,0.4)+
  draw_plot(p1a,0,0.35,1,0.65)+
  draw_label("A",0.025,0.95,size=11,fontface = "bold")+
  draw_label("B",0.025,0.4,size=11,fontface = "bold")
dev.off()

pdf("~/popvae/fig/ag1000g_inversion_scan_by_species.pdf",width=6.5,height=3.5,useDingbats = F)
ggdraw()+
  draw_plot(p2,0,0,1,0.4)+
  draw_plot(p1b,0,0.35,1,0.65)+
  draw_label("A",0.025,0.95,size=11,fontface = "bold")+
  draw_label("B",0.025,0.4,size=11,fontface = "bold")
dev.off()











