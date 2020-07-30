library(ggplot2);library(data.table);library(plyr);library(cowplot)
setwd("~/popvae/")

theme_set(theme_classic()+
            theme(axis.line = element_blank(),
                  #panel.grid.major = element_line(size=0.4,color="grey"),
                  axis.ticks = element_line(size=0.4),
                  strip.text = element_text(size=8),
                  legend.title=element_text(size=8,face="bold"),
                  legend.text=element_text(size=8),
                  axis.text=element_text(size=7),
                  strip.background = element_blank(),#rect(size=0.5,fill=NA),
                  axis.title=element_text(size=8)))

#LD decay
real <- fread("out/1kg/1kg_LD_decay_chr10:10000000.0-11000000.0_real.csv")
gen <- fread("out/1kg/1kg_LD_decay_chr10:10000000.0-11000000.0_gen.csv")
sim <- fread("out/1kg/1kg_LD_decay_chr10:10000000.0-11000000.0_sim.csv")

real$type <- "real"
gen$type <- "VAE"
sim$type <- "simulation"

df <- rbind(real,gen,sim)
df$LD <- round(df$LD,digits=3)
df <- subset(df,!is.na(LD))
df <- subset(df,dist<1e6)
df$bin <- cut(df$dist,breaks = 20)
binmeans <- ddply(df,.(bin,type),summarize,
                  mean=mean(LD),
                  median=median(LD),
                  high=quantile(LD,0.9),
                  low=quantile(LD,0.1),
                  x=median(dist))
#binmeans <- subset(binmeans,x!=max(binmeans$x[binmeans$type=="real"]))
#pdf("fig/1kg_LD_means.pdf",width=2.75,height=2.25)

LDplot <- ggplot(data=binmeans,aes(x=x/1e6,y=mean,col=type))+
  scale_color_manual(values=c("goldenrod","chartreuse4","steelblue3"),name="Genotypes",guide=F)+
  scale_fill_manual(values=c("goldenrod","chartreuse4","steelblue3"),name="Genotypes",guide=F)+
  ylab(expression(Linkage~Disequilibrium~(italic(R^2))))+
  xlab("Distance (Mbp)")+
  geom_line(lwd=0.35)+
  #geom_line(aes(y=high),linetype=2)+
  #geom_line(aes(y=low),linetype=2)+
  geom_point(size=1.2,shape=21,color="white",aes(fill=type),stroke=0.5)
  #ylim(0,0.1)
print(LDplot)
#dev.off() 


#comparing LD decay in 4 regions
bin1 <- "chr10:10000000.0-11000000.0"
bin2 <- "chr10:25000000.0-26000000.0"
bin3 <- "chr22:25000000.0-26000000.0"
bin4 <- "chr22:30000000.0-31000000.0"
for(i in c(bin1,bin2,bin3,bin4)){
  for(j in c("real","gen","sim")){
    f <- paste0("out/1kg/1kg_LD_decay_",i,"_",j,".csv")
    print(f)
    a <- fread(f)
    a$bin <- i
    a$type <- j
    a <- subset(a,!is.na(LD))
    if(i==bin1 & j=="real") b <- a else b <- rbind(b,a)
  }
}

b$window <- cut(b$dist,breaks = 25)

library(doMC)
registerDoMC(cores=12)
a <- ddply(b,.(bin,type,window),summarize,.parallel=TRUE,
           mean=mean(LD,na.rm=T),
           low=quantile(LD,0.75),
           high=quantile(LD,0.25),
           median=median(LD),
           x=median(dist)) 
a$type2 <- factor(a$type,levels=c("real","gen","sim"),labels=c("real","VAE decoder","coalescent\nsimulation"))
a$bin <- gsub("\\.0","",a$bin)

pdf("fig/1kg_LD_curve_window_comparisons.pdf",width=6.5,height=4,useDingbats = F)
ggplot(data=a,aes(x=x/1e6,y=median,col=type2,fill=type2))+
  facet_wrap(~bin,scales="free")+
  theme(#legend.position = c(0.875,0.825),
        axis.line.x=element_blank(),
        axis.ticks.length = unit(2,"mm"))+
  scale_color_manual(values=c("goldenrod","chartreuse4","steelblue3"),name="Genotypes")+
  scale_fill_manual(values=c("goldenrod","chartreuse4","steelblue3"),name="Genotypes")+
  ylab(expression(Linkage~Disequilibrium~(italic(R^2))))+
  xlab("Distance (Mbp)")+
  geom_line(lwd=0.4)+
  geom_hline(yintercept = 0,lwd=0.4)+
  geom_vline(xintercept = 0,lwd=0.4)+
  annotate("segment",x=0.15,xend=0.825,y=Inf,yend=Inf,color="black",lwd=0.7)+
  #geom_line(aes(y=high),linetype=2)+
  #geom_line(aes(y=low),linetype=2)+
  geom_point(size=1.25,shape=21,color="white",stroke=0.35)+
  guides(fill=guide_legend(keyheight = unit(4,"mm")),
         col=guide_legend(keyheight = unit(4,"mm")))
dev.off()

ggplot(data=subset(b,b$dist<1e5 & b$bin==b$bin[1]),aes(x=LD,fill=type))+
  facet_wrap(~bin)+
  xlim(0,0.5)+
  geom_density(alpha=0.5)



#PCA
a <- fread("out/1kg/1kg_decoder_PCA.csv")
b <- fread("out/1kg/1kg_sim_PCA.txt")
head(a)
head(b)
a <- arrange(a,pop)
b <- arrange(b,pop)
c <- data.frame(PC1=c(a$realPC1,a$genPC1,b$PC1),
                PC2=c(a$realPC2,a$genPC2,b$PC2),
                type=c(rep("real genotypes",150),rep("VAE decoder genotypes",150),rep("simulated genotypes",150)),
                pop=c(a$pop,a$pop,a$pop))

#pdf("fig/1kg_pca.pdf",width=6,height=2,useDingbats = F)
pcaplot <- ggplot(data=c,aes(x=PC1,y=PC2,fill=pop))+
  theme(axis.line = element_blank(),
        panel.border = element_rect(color="black",fill=NA))+
  facet_wrap(~type)+
  scale_fill_brewer(palette = "RdYlBu",name="Population")+
  geom_vline(xintercept=0,lwd=0.25)+geom_hline(yintercept=0,lwd=0.25)+
  geom_point(shape=21,stroke=0.3)+
  guides(fill=guide_legend(override.aes = list(size=4),keyheight = unit(4,"mm")))
print(pcaplot)
#dev.off()

#SFS
sfs <- fread("out/1kg/1kg_sfs.csv")
sfs <- subset(sfs,bin>2 & bin<100)
sfs$bin <- sfs$bin/100
sfs <- reshape::melt(sfs,id.vars="bin")
sfs$variable <- factor(sfs$variable,levels=c("real","simulation","VAE"))
#sfs$value[sfs$variable=="simulation"] <- sfs$value[sfs$variable=="simulation"]*(117781/48350) #simulation SFS seems to be off by a factor of ~2.5(?)

#pdf("fig/1kg_sfs.pdf",width=3,height=2.5,useDingbats = F)
sfsplot <- ggplot(data=sfs,aes(x=bin,y=value,color=variable))+
  scale_y_log10()+
  theme(legend.box.margin=margin(0,0,0,-15))+
  scale_color_manual(values=c("goldenrod","chartreuse4","steelblue3"),name="Genotypes")+
  xlab("Allele Frequency")+ylab("Sites")+
  annotation_logticks(sides="l",long=unit(2.25,"mm"),mid=unit(1.75,"mm"),short=unit(1,"mm"))+
  geom_line()+
  guides(color=guide_legend(keyheight = unit(4,"mm")))
print(sfsplot)
#dev.off()


pdf("fig/1kg_pca+sfs+ld.pdf",width=6,height=3.5,useDingbats = F)
ggdraw()+
  draw_plot(pcaplot,0,0.5,1,0.5)+
  draw_plot(LDplot,0.05,0,0.4,0.5)+
  draw_plot(sfsplot,0.455,0,0.55,0.5)+
  draw_label("A",0.025,0.95)+
  draw_label("B",0.025,0.48)+
  draw_label("C",0.465,0.48)
dev.off()




