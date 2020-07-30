library(ggplot2);library(reshape);library(plyr);
library(data.table);library(raster);library(broom);library(cowplot)

setwd("~/popvae_dev/")
theme_set(theme_classic()+theme(axis.title=element_text(size=7),
                                axis.text=element_text(size=6),
                                legend.text=element_text(size=7),
                                legend.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                strip.background = element_blank()))

#function for plotting the first eight PC axes as a grid (pc should be a data frame with columns 'PC1 ... PC8 ... <fill>').
pal <- c(RColorBrewer::brewer.pal(12,"Paired"),"black","grey","violet","navyblue")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette


########################################
################# HGDP #################
########################################
a <- fread("out/hgdp/hgdp_latent_coords.txt",header=T)
names(a) <- c("LD1","LD2","sampleID")
pc <- fread("out/hgdp/hgdp_defaults_pca.txt",data.table = F)
pc <- pc[,c(1,2,ncol(pc))]
names(pc) <- c("LD1","LD2","sampleID")
b <- fread("data/hgdp/hgdp_sample_data.txt")
b$region <- factor(b$region,labels=c("Africa","Americas","Central/South Asia","East Asia","Europe","Middle East","Oceania"),
                            levels=c("AFRICA","AMERICA","CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST","OCEANIA"))
c <- merge(a,b,by="sampleID")
d <- merge(b,pc,by="sampleID")
c$method <- "VAE"
d$method <- "PCA"
d$LD1 <- -d$LD1
d$LD2 <- d$LD2
e <- rbind(c,d)

labels <- ddply(e,.(population,method),summarize,ld1=mean(LD1),ld2=mean(LD2))

pc <- fread("out/hgdp/hgdp_defaults_pca.txt")
pc <- merge(pc,b,by="sampleID")
plot_multipanel_pc <- function(pc,fill,pal,nrow=2){
  p1 <- ggplot(data=pc,aes_string(x="PC2",y="-PC1",fill=fill))+
    theme(axis.title.y=element_text(margin=margin(0,-2,0,0)),
          plot.background = element_blank())+
    geom_point(shape=21,stroke=0.05,size=1.75)+
    coord_cartesian(clip="off")+
    scale_fill_manual(values = pal,guide=F)
  p2 <- ggplot(data=pc,aes_string(x="PC3",y="PC4",fill=fill))+
    theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
    geom_point(shape=21,stroke=0.05,size=1.75)+
    coord_cartesian(clip="off")+
    scale_fill_manual(values = pal,guide=F)
  p3 <- ggplot(data=pc,aes_string(x="PC5",y="PC6",fill=fill))+
    theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
    geom_point(shape=21,stroke=0.05,size=1.75)+
    coord_cartesian(clip="off")+
    scale_fill_manual(values = pal,guide=F)
  p4 <- ggplot(data=pc,aes_string(x="PC7",y="PC8",fill=fill))+
    theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
    geom_point(shape=21,stroke=0.05,size=1.75)+
    coord_cartesian(clip="off")+
    scale_fill_manual(values = pal,guide=F)
  #p <- plot_grid(p1,p2,p3,p4,nrow = nrow)
  p <- ggdraw()+
    draw_plot(p2,0,0,0.34,0.435)+
    draw_plot(p3,0.33,0,0.34,0.435)+
    draw_plot(p4,0.66,0,0.34,0.435)+
    draw_plot(p1,0,0.4,1,0.6)
  return(p)
}
p1 <- plot_multipanel_pc(pc,"region",cbPalette)
p2 <- ggplot(c,aes(x=-LD2,y=LD1,fill=region))+
  theme(legend.spacing.y = unit(0,"mm"),
        legend.position="bottom",
        legend.box.margin=margin(-10,0,0,0))+
  geom_point(stroke=0.05,color="black",alpha=0.8,shape=21,size=1.75)+
  scale_fill_manual(values=cbPalette,name="Region")+
  guides(fill=guide_legend(override.aes = list(size=4,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "top",
                           ncol=3))

ggplot(c,aes(x=-LD2,y=LD1,fill=longitude))+
  theme(legend.position=c(0.07,0.8),
        legend.background = element_blank())+
  geom_point(stroke=0.1,color="black",alpha=0.8,shape=21,size=1.75)+
  scale_fill_distiller(palette = "YlGnBu",name="Longitude")


pdf("fig/hgdp/hgdp_pca_v_vae5.pdf",width=6.5,height=3.5,useDingbats = F)
ggdraw()+
  draw_plot(p1,0,0,0.5,0.95)+
  draw_plot(p2,0.5,0,0.5,0.95)+
  draw_text("PCA",size=10,0.25,0.95)+
  draw_text("VAE",size=10,0.75,0.95)
dev.off()

# tall version
fill="region"
pc1 <- ggplot(data=pc,aes_string(x="PC2",y="-PC1",fill=fill))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)),
        plot.background = element_blank())+
  ggtitle("PCA")+
  geom_point(shape=21,stroke=0.05,size=1.75)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)
pc2 <- ggplot(data=pc,aes_string(x="PC3",y="PC4",fill=fill))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(shape=21,stroke=0.05,size=1.75)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)
pc3 <- ggplot(data=pc,aes_string(x="PC5",y="PC6",fill=fill))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(shape=21,stroke=0.05,size=1.75)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)
pc4 <- ggplot(data=pc,aes_string(x="PC7",y="PC8",fill=fill))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(shape=21,stroke=0.05,size=1.75)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)
pc5 <- ggplot(data=pc,aes_string(x="PC9",y="PC10",fill=fill))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(shape=21,stroke=0.05,size=1.75)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)
p2 <- ggplot(c,aes(x=LD2,y=-LD1,fill=region))+
  theme(legend.spacing.y = unit(0,"mm"),
        legend.box.margin=margin(-10,0,0,0))+
  ggtitle("VAE")+
  geom_point(stroke=0.05,color="black",alpha=0.8,shape=21,size=1.75)+
  scale_fill_manual(values=cbPalette,name="Region")+
  guides(fill=guide_legend(override.aes = list(size=4,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "top",
                           ncol=1))
  


pdf("fig/hgdp/hgdp_pca_v_vae_tall.pdf",width=6.5,height=7,useDingbats = F)
ggdraw()+
  draw_plot(pc1,0,0.5,0.66,0.5)+
  draw_plot(pc2,0.66,0.5+2*(.5*.33),0.33,0.5*.33)+
  draw_plot(pc3,0.66,0.5+(.5*.33),0.33,0.5*.33)+
  draw_plot(pc4,0.66,0.5,0.33,0.5*.33)+
  draw_plot(p2,0,0,1,0.5)
dev.off()


#large plot with population labels
library(ggrepel)
pop_centroids <- ddply(c,.(population,region),summarize,LD1=mean(LD1),LD2=mean(LD2))
pdf("fig/hgdp/hgdp_vae_defaults_large_poplabels.pdf",width=6,height=5,useDingbats=F)
ggplot(c,aes(x=-LD2,y=LD1,fill=region))+
  theme(legend.spacing.y = unit(0,"mm"),
        legend.position="bottom",
        legend.box.margin=margin(-10,0,0,0))+
  geom_point(stroke=0.05,color="black",alpha=0.8,shape=21,size=1.75)+
  geom_text_repel(data=pop_centroids,aes(label=population),size=2.5)+
  scale_fill_manual(values=cbPalette,name="Region")+
  guides(fill=guide_legend(override.aes = list(size=4,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "top",
                           ncol=3))
dev.off()

# #tree of pairwise VAE distances
# library(ape);library(sp)
# dists <- spDists(x=as.matrix(c[,c("LD1","LD2")]))
# rownames(dists) <- c$region
# tree <- nj(dists)
# tree <- root(tree,"Africa")
# plot(tree,cex=0.3)


#summarized one latent dimension a map
load("~/locator/locator_hgdp/cntrymap.Rdata")
map <- crop(map,c(-170,170,-50,71))
a <- fread("out/hgdp/hgdp_6x128_1LD_1e5snps_latent_coords.txt",header=T)
names(a) <- c("LD1","sampleID")
a <- merge(a,b,by="sampleID")
d <- ddply(a,.(longitude,latitude,population),summarize,n=length(LD1),mean_LD1=mean(LD1))
p <- ggplot()+coord_map("mollweide")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        #axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.15,color="grey60"),
        legend.position=c(0.075,0.55),
        #legend.background = element_blank(),
        legend.spacing = unit(0,"mm"))+
  scale_fill_distiller(palette = "YlGnBu",name="Mean LD1")+
  scale_size_continuous(name="Samples",breaks=c(10,20,30,40))+
  geom_polygon(data=map,aes(x=long,y=lat,group=group),lwd=0.1,col="white",fill="grey60")+
  geom_point(data=d,aes(x=longitude,y=latitude,fill=mean_LD1,group=population,size=n),
             shape=21,stroke=0.2,alpha=0.8)+ #can flip LD1 v LD2 here
  guides(fill=guide_colorbar(barheight = unit(18,"mm"),barwidth = unit(3.5,"mm")),
         size=guide_legend(keyheight = unit(2,"mm"),override.aes = list(stroke=0.5)))


regionmeans <- ddply(a,.(region),summarize,LD1=mean(LD1))

dplot <- ggplot(data=a,aes(x=LD1,group=region,fill=region))+
  theme(legend.position = "none")+
  ylab("Density")+xlab("LD1")+
  scale_fill_manual(values=cbPalette,name="Region")+
  geom_text_repel(data=regionmeans,aes(x=LD1,group=region,y=0,label=region),nudge_y = 5,nudge_x=-0.175,size=2.5)+
  geom_density(lwd=0.25,alpha=1)
  

pdf("fig/hgdp/hgdp_LD1_map.pdf",useDingbats = F,width=5.8,height=3.5)
plot_grid(dplot,p,rel_heights=c(0.35,0.7),ncol=1)
dev.off()

ggplot(data=a,aes(x=LD1,fill=region))+
  scale_fill_manual(values=pal,name="Region")+
  geom_density()

###as above but with PC1
load("~/locator/locator_hgdp/cntrymap.Rdata")
map <- crop(map,c(-170,170,-50,71))
d <- ddply(pc,.(longitude,latitude,population),summarize,n=length(PC1),mean_PC1=mean(PC1))
p <- ggplot()+coord_map("mollweide")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        #axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.15,color="grey60"),
        legend.position=c(0.075,0.55),
        #legend.background = element_blank(),
        legend.spacing = unit(0,"mm"))+
  scale_fill_distiller(palette = "YlGnBu",name="Mean PC1")+
  scale_size_continuous(name="Samples",breaks=c(10,20,30,40))+
  geom_polygon(data=map,aes(x=long,y=lat,group=group),lwd=0.1,col="white",fill="grey60")+
  geom_point(data=d,aes(x=longitude,y=latitude,fill=mean_PC1,group=population,size=n),
             shape=21,stroke=0.2,alpha=0.8)+ #can flip LD1 v LD2 here
  guides(fill=guide_colorbar(barheight = unit(18,"mm"),barwidth = unit(3.5,"mm")),
         size=guide_legend(keyheight = unit(2,"mm"),override.aes = list(stroke=0.5)))

regionmeans <- ddply(pc,.(region),summarize,PC1=mean(PC1))

dplot <- ggplot(data=pc,aes(x=PC1,group=region,fill=region))+
  theme(legend.position = "none")+
  ylab("Density")+xlab("LD1")+
  scale_fill_manual(values=cbPalette,name="Region")+
  ggrepel::geom_text_repel(data=regionmeans,aes(x=PC1,group=region,y=0,label=region),nudge_y = 5,nudge_x=-0.175,size=2.5)+
  geom_density(lwd=0.25,alpha=1)


pdf("fig/hgdp/hgdp_PC1_map.pdf",useDingbats = F,width=5.8,height=3.5)
plot_grid(dplot,p,rel_heights=c(0.35,0.7),ncol=1)
dev.off()

ggplot(data=a,aes(x=LD1,fill=region))+
  scale_fill_manual(values=pal,name="Region")+
  geom_density()


#animated gif of (centered) predictions during training
library(gganimate)
a2 <- fread("out/hgdp/hgdp_s94019_p20_training_preds.txt")
names(a2) <- c("LD1","LD2","sampleID","epoch")
b2 <- fread("data/hgdp/hgdp_sample_data.txt")
c2 <- merge(a2,b2,by="sampleID")
c2 <- ddply(c2,.(epoch),function(k) { #center predictions for prettier plots
  k$LD1 <- k$LD1-mean(k$LD1)
  k$LD2 <- k$LD2-mean(k$LD2)
  return(k)
})
ranges <- ddply(c2,.(epoch),summarize,r2=max(LD2)-min(LD2),r1=max(LD1)-min(LD1))
dropranges <- subset(ranges,r2==max(ranges$r2)|r1==max(ranges$r1))
c2 <- subset(c2,!(epoch %in% dropranges$epoch))
c2$region <- factor(c2$region,labels=c("Africa","America","Central/South Asia","East Asia","Europe","Middle East","Oceania"),
                   levels=c("AFRICA","AMERICA","CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST","OCEANIA"))

ggplot(c2,aes(x=-LD2,y=LD1,fill=region))+
  geom_point(shape=21,stroke=0.2,size=3)+
  scale_fill_manual(values = cbPalette)+
  transition_manual(epoch)+
  guides(fill=guide_legend(override.aes = list(size=4)))+
  labs(title='Epoch: {frame}')
anim_save("out/hgdp_s94019_p20_training_preds.gif")

c3 <- subset(c2,epoch %in% c(0,5,10,20,30,40,60,80,100,160,320,500))
pdf("fig/hgdp_training_preds_by_epoch.pdf",useDingbats = F,width=6.5,height=5)
ggplot(c3,aes(x=-LD2,y=-LD1,fill=region))+
  theme(legend.position = c(0.88,0.15))+
  facet_wrap(~epoch,scales = "free",labeller = function(e) lapply(e,function(i)paste0("epoch ",i)))+
  geom_point(shape=21,stroke=0.2)+
  scale_fill_manual(values = cbPalette,name="Region")+
  guides(fill=guide_legend(override.aes = list(size=4),keyheight = unit(0,"mm")))
dev.off()


##########################################
################# ag1000g ################
##########################################
a <- fread("out/ag1000g/ag1000g_defaults_latent_coords.txt",header=T)
names(a) <- c("LD1","LD2","sampleID")
pc <- fread("out/ag1000g/ag1000g_defaults_pca.txt")[,c(1,2,21)]
names(pc) <- c("LD1","LD2","sampleID")
b <- fread("~/locator/data/ag1000g/anopheles_samples_sp.txt")
b$species[is.na(b$species)] <- "Uncertain"
#b$sampleID <- b$ox_code
c <- merge(a,b,by="sampleID")
d <- merge(b,pc,by="sampleID")
c$method <- "VAE"
d$method <- "PCA"
e <- rbind(c,d)

pc <- fread("out/ag1000g/ag1000g_512x4_pca.txt")
pc <- merge(pc,b,"sampleID")
pc$country[pc$country=="France"] <- "Mayotte"
pc$country[pc$country=="Gambia, The"] <- "The Gambia"

p1 <- ggplot(data=pc,aes(x=PC1,y=PC2,fill=country,shape=species))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(stroke=0.05,size=2)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)+
  scale_shape_manual(values=c(21,22,23),guide=F)
p2 <- ggplot(data=pc,aes(x=PC3,y=PC4,fill=country,shape=species))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(stroke=0.05,size=2)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)+
  scale_shape_manual(values=c(21,22,23),guide=F)
p3 <- ggplot(data=pc,aes(x=PC5,y=PC6,fill=country,shape=species))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(stroke=0.05,size=2)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)+
  scale_shape_manual(values=c(21,22,23),guide=F)
p4 <- ggplot(data=pc,aes(x=PC7,y=PC8,fill=country,shape=species))+
  theme(axis.title.y=element_text(margin=margin(0,-2,0,0)))+
  geom_point(stroke=0.05,size=2)+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = pal,guide=F)+
  scale_shape_manual(values=c(21,22,23),guide=F)
p <- ggdraw()+
  draw_plot(p1,0,0.4,1,0.6)+
  draw_plot(p2,0,0,0.35,0.43)+
  draw_plot(p3,0.33,0,0.35,0.43)+
  draw_plot(p4,0.66,0,0.35,0.43)


c$country[c$country=="France"] <- "Mayotte"
c$country[c$country=="Gambia, The"] <- "The Gambia"

p2 <- ggplot(c,aes(x=LD1,y=LD2,fill=country,shape=species))+
  theme(legend.spacing.y = unit(0,"mm"))+
  geom_point(stroke=0.1,color="black",alpha=0.8,size=2)+
  scale_fill_manual(values=pal,name="Country")+
  scale_shape_manual(values=c(21,22,23),name="Species",labels=c(expression(italic(A.~coluzzii)),
                                                                expression(italic(A.~gambiae)),
                                                                "unknown"))+
  guides(fill=guide_legend(override.aes = list(size=4,shape=21),
                           keyheight = unit(3,"mm"),
                           title.position = "top",
                           ncol=1),
         shape=guide_legend(override.aes = list(size=4,stroke=0.8),
                            keyheight = unit(1,"mm"),
                            ncol=1))

pdf("fig/ag1000g/ag1000g_pca_v_vae.pdf",width=6.5,height=3.25,useDingbats = F)
ggdraw()+
  draw_plot(p,0,0,0.425,0.95)+
  draw_plot(p2,0.425,0,0.575,0.95)+
  draw_text("PCA",size=10,0.25,0.95)+
  draw_text("VAE",size=10,0.65,0.95)
dev.off()

#animated gif of (centered) predictions during training
library(gganimate)
a2 <- fread("out/ag1000g_phase2_3R_1e5snps_2x128_training_preds.txt")
names(a2) <- c("LD1","LD2","sampleID","epoch")
b2 <- fread("data/ag1000g/anopheles_samples_sp.txt")
c2 <- merge(a2,b2,by="sampleID")
c2 <- ddply(c2,.(epoch),function(k) { #center predictions for prettier plots
  k$LD1 <- k$LD1-mean(k$LD1)
  k$LD2 <- k$LD2-mean(k$LD2)
  return(k)
})
ranges <- ddply(c2,.(epoch),summarize,r2=max(LD2)-min(LD2),r1=max(LD1)-min(LD1))
dropranges <- subset(ranges,r2==max(ranges$r2)|r1==max(ranges$r1))
c2 <- subset(c2,!(epoch %in% dropranges$epoch))
ggplot(c2,aes(x=LD1,y=LD2,fill=population))+
  geom_point(shape=21,stroke=0.2,size=3)+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(12,"Paired"),"black","grey","violet","navyblue"),name="population")+
  transition_manual(epoch)+
  labs(title='Epoch: {frame}')
anim_save("out/ag1000g.gif")

