library(ggplot2);library(reshape);library(plyr);library(tsne)
library(data.table);library(raster);library(broom);library(cowplot)

setwd("~/popvae_dev/")
theme_set(theme_classic()+theme(axis.title=element_text(size=7),
                                axis.text=element_text(size=6),
                                legend.text=element_text(size=7),
                                legend.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                title = element_text(size=8),
                                strip.background = element_blank()))

pal <- c(RColorBrewer::brewer.pal(12,"Paired"),"black","grey","violet","navyblue")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette

#load run times for HGDP
rt_1e5 <- fread("out/hgdp/hgdp_defaults_1e5snps_runtimes.txt")
rt_1e5$nsnps <- 1e5
rt_1e5$method <- c("VAE","PCA")

rt_1e4 <- fread("out/hgdp/hgdp_defaults_1e4snps_runtimes.txt")
rt_1e4$nsnps <- 1e4
rt_1e4$method <- c("VAE","PCA")

#get umap and t-SNE run times
library(umap)
pc1 <- fread("out/hgdp/hgdp_defaults_1e5snps_pca.txt")
t1 <- Sys.time()
um1 <- umap(pc1[,1:15])
diff1 <- Sys.time()-t1

t1 <- Sys.time()
tsne1 <- tsne(pc1[,1:15])
diff2 <- Sys.time()-t1

pc2 <- fread("out/hgdp/hgdp_defaults_1e4snps_pca.txt")
t1 <- Sys.time()
um2 <- umap(pc2[,1:15])
diff3 <- Sys.time()-t1

t1 <- Sys.time()
tsne2 <- tsne(pc2[,1:15])
diff4 <- Sys.time()-t1

#plot for a reality check
umap1 <- as.data.frame(um1$layout)
umap2 <- as.data.frame(um2$layout)
tsne1 <- as.data.frame(tsne1)
tsne2 <- as.data.frame(tsne2)
umap1$method <- "UMAP"
umap2$method <- "UMAP"
tsne1$method <- "t-SNE"
tsne2$method <- "t-SNE"
umap1$snps <- "100,000 SNPs"
umap2$snps <- "10,000 SNPs"
tsne1$snps <- "100,000 SNPs"
tsne2$snps <- "10,000 SNPs"
umap1$sampleID <- pc1$sampleID
umap2$sampleID <- pc2$sampleID
tsne1$sampleID <- pc1$sampleID
tsne2$sampleID <- pc2$sampleID
names(umap1) <- names(tsne1)
names(umap2) <- names(tsne2)
pd <- rbind(umap1,umap2,tsne1,tsne2)
sd <- fread("data/hgdp/hgdp_sample_data.txt")
pd <- merge(pd,sd,by="sampleID")
pd$region <- factor(pd$region,labels=c("Africa","America","Central/South Asia","East Asia","Europe","Middle East","Oceania"),
                   levels=c("AFRICA","AMERICA","CENTRAL_SOUTH_ASIA","EAST_ASIA","EUROPE","MIDDLE_EAST","OCEANIA"))

pdf("fig/hgdp_umap_tsne_plots_15pc.pdf",useDingbats = F,width=6,height=3.5)
ggplot(data=pd,aes(x=V1,y=V2,fill=region))+
  facet_wrap(method~snps,scales="free")+
  theme(axis.title = element_blank())+
  geom_point(shape=21,stroke=0.05)+
  scale_fill_manual(values = cbPalette,name="Region")+
  guides(fill=guide_legend(keyheight = unit(0,"mm"),override.aes = list(size=3)))
dev.off()

print(diff1)
print(diff2)

#get normalized pairwise spatial and latent space distances for comparisons
vae <- fread("out/hgdp/hgdp_defaults_latent_coords.txt",header = T)
names(vae)[1:2] <- c("V1","V2")
vae$method <- "VAE"
vae$snps <- "100,000 SNPs"
pca <- fread("out/hgdp/hgdp_defaults_1e5snps_pca.txt",header=T)
pca <- pca[,c("PC1","PC2","sampleID")]
names(pca)[1:2] <- c("V1","V2")
pca$snps <- "100,000 SNPs"
pca$method <- "PCA"

df <- rbind(vae,umap1,tsne1,pca)
df <- merge(df,sd,by="sampleID")
df <- subset(df,region %in% c("CENTRAL_SOUTH_ASIA",'EUROPE','MIDDLE_EAST','EAST_ASIA')) #looking only at samples from Eurasia where pairwise distance is at least somewhat informative

dist_pca <- c(as.matrix(dist(df[df$method=="PCA",c("V1","V2")])))
dist_pca <- dist_pca/max(dist_pca)
dist_vae <- c(as.matrix(dist(df[df$method=="VAE",c("V1","V2")])))
dist_vae <- dist_vae/max(dist_vae)
dist_umap <- c(as.matrix(dist(df[df$method=="UMAP",c("V1","V2")])))
dist_umap <- dist_umap/max(dist_umap)
dist_tsne <- c(as.matrix(dist(df[df$method=="t-SNE",c("V1","V2")])))
dist_tsne <- dist_tsne/max(dist_tsne)

#reconciling population names
# ol_dists <- fread("data/hgdp/Ramachandran_overland_distances.txt")
# ol_dists$pop1[ol_dists$pop1=="Italian"] <- "BergamoItalian"
# ol_dists$pop1[ol_dists$pop1=="BantuSouthWest"] <- "BantuSouthAfrica"
# ol_dists$pop1[ol_dists$pop1=="BantuSouthEast"] <- "BantuKenya"
# ol_dists$pop1[ol_dists$pop1=="BiakaPygmy"] <- "Biaka"
# ol_dists$pop1[ol_dists$pop1=="MbutiPygmy"] <- "Mbuti"
# ol_dists$pop1[ol_dists$pop1=="Han-NChina"] <- "NorthernHan"
# ol_dists$pop1[ol_dists$pop1=="Mongola"] <- "Mongolian"
# ol_dists$pop1[ol_dists$pop1=="Melanesian"] <- "Bougainville"
# ol_dists$pop2[ol_dists$pop2=="Italian"] <- "BergamoItalian"
# ol_dists$pop2[ol_dists$pop2=="BantuSouthWest"] <- "BantuSouthAfrica"
# ol_dists$pop2[ol_dists$pop2=="BantuSouthEast"] <- "BantuKenya"
# ol_dists$pop2[ol_dists$pop2=="BiakaPygmy"] <- "Biaka"
# ol_dists$pop2[ol_dists$pop2=="MbutiPygmy"] <- "Mbuti"
# ol_dists$pop2[ol_dists$pop2=="Han-NChina"] <- "NorthernHan"
# ol_dists$pop2[ol_dists$pop2=="Mongola"] <- "Mongolian"
# ol_dists$pop2[ol_dists$pop2=="Melanesian"] <- "Bougainville"
# ol_dists$pair <- paste(ol_dists$pop1,ol_dists$pop2)
# df$population[df$population %in% c("PapuanSepik","PapuanHighlands")] <- "Papuan"
# df$population[df$population %in% c("Surui")] <- "Karitiana"
# 
# pops <- df[df$method=="VAE"]$population
# pop_pairs <- matrix(nrow=length(pops),ncol=length(pops))
# for(i in 1:length(pops)){
#   for(j in 1:length(pops)){
#     pop_pairs[i,j] <- paste(pops[i],pops[j])
#   }
# }
# pop_pairs <- c(pop_pairs)
# pops2 <- unique(c(ol_dists$pop1,ol_dists$pop2))
# pops2[!pops2 %in% unique(pops)]
# unique(df$population)[!unique(df$population) %in% pops2]
# 
# ol_dists_sorted <- rep(1.0,length(pop_pairs))
# for(i in 1:length(pop_pairs)){
#   ol_dists_sorted[i] <- ol_dists$corrected[ol_dists$pair==pop_pairs[i]][1]
#   print(i)
# }
# dist_geo <- ol_dists_sorted

dist_geo <- c(spDists(as.matrix(df[df$method=="VAE",c("longitude","latitude")]),longlat = T))
dist_geo <- dist_geo/max(dist_geo)

cor_pca <- cor(dist_pca,dist_geo)^2
cor_vae <- cor(dist_vae,dist_geo)^2
cor_umap <- cor(dist_umap,dist_geo)^2
cor_tsne <- cor(dist_tsne,dist_geo)^2
r2 <- c(cor_pca,cor_vae,cor_umap,cor_tsne)
r2 <- sapply(r2,function(e) paste0("italic(R)^2==",round(e,3)))

pd <- data.frame(geographic=dist_geo,
                 UMAP=dist_umap,
                 tSNE=dist_tsne,
                 VAE=dist_vae,
                 PCA=dist_pca)
mpd <- reshape::melt(pd,id.vars="geographic")
ld <- data.frame(x=c(0,1),y=c(0,1))
labels <- data.frame(x=c(0.18,0.18,0.78,0.78),
                     y=rep(1.04,4),
                     r2=r2,
                     variable=c("PCA","VAE","UMAP","tSNE"))

pdf("fig/dist_cor_comparisons_1row.pdf",width=6.5,height=2,useDingbats = F)
p <- ggplot(data=mpd,aes(x=geographic,y=value))+
  coord_equal()+
  theme(panel.border = element_rect(fill=NA,size=1),axis.line = element_blank())+
  #ggtitle("Latent Space and Geographic Distance, by Method\nHuman Genotypes from Europe and Asia")+
  facet_wrap(~variable,nrow=1)+
  xlab("Geographic Distance")+ylab("Latent Space Distance")+
  scale_fill_distiller(palette = "YlGnBu",name="Count")+
  ylim(0,1.06)+
  stat_bin_hex()+
  geom_line(data=ld,aes(x=x,y=y))+
  geom_text(data=labels,aes(x=x,y=y,label=r2),parse=T,size=2.5)+
  guides(fill=guide_colorbar(barwidth = unit(3,"mm"),barheight = unit(18,"mm")))
print(p)
dev.off()




##################### ag1000g ########################
setwd("~/popvae_dev/")
theme_set(theme_classic()+theme(axis.title=element_text(size=7),
                                axis.text=element_text(size=6),
                                legend.text=element_text(size=7),
                                legend.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                strip.background = element_blank()))

pal <- c(RColorBrewer::brewer.pal(12,"Paired"),"black","grey","violet","navyblue")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette

#load run times for HGDP
rt_1e5 <- fread("out/ag1000g/ag1000g_1e5snps_runtimes.txt")
rt_1e5$nsnps <- 1e5
rt_1e5$method <- c("VAE","PCA")

rt_1e4 <- fread("out/ag1000g/ag1000g_1e4snps_runtimes.txt")
rt_1e4$nsnps <- 1e4
rt_1e4$method <- c("VAE","PCA")

#get umap and t-SNE run times
library(umap)
pc1 <- fread("out/ag1000g/ag1000g_1e5snps_pca.txt")
t1 <- Sys.time()
#settings <- umap.defaults
#settings$init <- pc1[,1:2]
um1 <- umap(pc1[,1:15])
diff1 <- Sys.time()-t1

t1 <- Sys.time()
tsne1 <- tsne(pc1[,1:15])
diff2 <- Sys.time()-t1

pc2 <- fread("out/ag1000g/ag1000g_1e4snps_pca.txt")
t1 <- Sys.time()
um2 <- umap(pc2[,1:15])
diff3 <- Sys.time()-t1

t1 <- Sys.time()
tsne2 <- tsne(pc2[,1:15])
diff4 <- Sys.time()-t1

#plot for a reality check
umap1 <- as.data.frame(um1$layout)
umap2 <- as.data.frame(um2$layout)
tsne1 <- as.data.frame(tsne1)
tsne2 <- as.data.frame(tsne2)
umap1$method <- "UMAP"
umap2$method <- "UMAP"
tsne1$method <- "t-SNE"
tsne2$method <- "t-SNE"
umap1$snps <- "100,000 SNPs"
umap2$snps <- "10,000 SNPs"
tsne1$snps <- "100,000 SNPs"
tsne2$snps <- "10,000 SNPs"
umap1$sampleID <- pc1$sampleID
umap2$sampleID <- pc2$sampleID
tsne1$sampleID <- pc1$sampleID
tsne2$sampleID <- pc2$sampleID
#names(umap1)[1:2] <- c("X1","X2")
#names(umap2)[1:2] <- c("X1","X2")
pd <- rbind(umap1,umap2,tsne1,tsne2)
sd <- fread("data/ag1000g/anopheles_samples_sp.txt")
pd <- merge(pd,sd,by="sampleID")
pd$species[is.na(pd$species)] <- "unknown"
pd$country[pd$country=="France"] <- "Mayotte"
pd$country[pd$country=="Gambia, The"] <- "The Gambia"
pd$species[pd$species=="coluzzi"] <- "coluzzii"


pdf("fig/ag1000g_umap_tsne_plots_15pcs.pdf",useDingbats = F,width=6,height=3.5)
ggplot(data=pd,aes(x=V1,y=V2,fill=country,shape=species))+
  facet_wrap(method~snps,scales="free")+
  theme(legend.spacing.y = unit(0,"mm"))+
  geom_point(stroke=0.1,color="black",alpha=0.8,size=2)+
  scale_fill_manual(values=pal,name="Country")+
  scale_shape_manual(values=c(21,22,23),name="Species",labels=c(expression(italic(A.~coluzzi)),
                                                                expression(italic(A.~gambiae)),
                                                                "unknown"))+
  guides(fill=guide_legend(override.aes = list(size=4,shape=21),
                           keyheight = unit(1,"mm"),
                           title.position = "top",
                           ncol=1),
         shape=guide_legend(override.aes = list(size=4,stroke=0.8),
                            keyheight = unit(1,"mm"),
                            ncol=1))
dev.off()

