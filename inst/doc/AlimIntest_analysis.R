## ----library_loading, include=FALSE--------------------------------------
user=Sys.info()[["user"]]
require(rmarkdown)
require(markdown)
library(ade4)
library(RColorBrewer)
library(BiotypeR)
library(gridExtra)
library(reshape2)
library(vegan)
library(knitr)
library(ggplot2)
library(scales)
library(dplyr)
#source("src/mclapply.hack.R")

## ----knit_options, include=FALSE-----------------------------------------
opts_chunk$set(fig.width=12, fig.height=8, warning=FALSE, 
message=FALSE, echo=FALSE, dev=c('png', 'pdf', 'win.metafile', 'tiff'), 
fig.cap=NULL, cache=TRUE, dpi=300)

## ----replicate-----------------------------------------------------------

library(AlimIntest)
data(sanger_454)
data(replicat_454)

sanger_454_meta             = unique(sanger_454[1:3])[order(unique(sanger_454[1:3])$samples),]
rownames(sanger_454_meta)   = 1:dim(sanger_454_meta)[1]
replicat_454_meta           = unique(replicat_454[1:2])[order(unique(replicat_454[1:2])$samples),]
rownames(replicat_454_meta) = 1:dim(replicat_454_meta)[1]


sanger_454$tax[which(sanger_454$conf < 0.8)]     <- "Unassigned"
replicat_454$tax[which(replicat_454$conf < 0.8)] <- "Unassigned"


replicat_454.freq = replicat_454 %>% 
                    group_by(samples,manufacturer,tax) %>% 
					summarise (n = n()) %>% 
					mutate(freq = n / sum(n))
					
					
sanger_454.freq   = sanger_454 %>% 
                    group_by(samples,technology,run,tax) %>% 
					summarise (n = n()) %>% 
					mutate(freq = n / sum(n))
					
sanger_454.cast   = dcast(sanger_454.freq, tax ~ samples, fill = 0)
sanger_454.melt   = melt(sanger_454.cast, id.vars=c("tax","A"))

replicat_454.cast = dcast(replicat_454.freq, tax + samples ~ manufacturer, fill=0 )



## ----replicate_analysis_plot, fig.cap="Technical replicates comparisons based on genera relative log10 proportion", fig.height=3.5, fig.width=11.5----

sanger_454.plot = ggplot(sanger_454.melt, aes(x=A+10^-5, y=value+10^-5, group=variable)) + 
                  geom_point(aes(color=variable), size=2, alpha=0.8) + 
                  scale_x_log10("Sanger method",limits=c(10^-5,1)) + 
                  scale_y_log10("454 GS FLX Ti method",limits=c(10^-5,1)) + 
                  theme_bw() + geom_abline(slope=1, linetype = 2) + 
				  scale_color_discrete("454 method run",labels=sanger_454_meta$run[-1] )



replicat_454.plot = ggplot(replicat_454.cast, aes(x=A, y=B)) + 
                    geom_point(aes(color=samples), size=2, alpha=0.8) +
					scale_x_log10("454 GS FLX Ti Manufacturer A",limits=c(10^-5,1)) + 
                    scale_y_log10("454 GS FLX Ti Manufacturer B",limits=c(10^-5,1)) + 
                    theme_bw() + geom_abline(slope=1, linetype = 2)

grid.arrange(replicat_454.plot, sanger_454.plot, ncol=2)


## ----alimintest_study_design, results="asis", comment=NA, echo=FALSE-----


data(alimintestData)
kable(format(head(alimintestData$metadata[1:6])))


## ----alimintest_scfa, results="asis", comment=NA, echo=FALSE-------------

kable(format(head(alimintestData$metadata[7:14])))


## ----alimintest_qPCR, results="asis", comment=NA, echo=FALSE-------------

kable(format(head(alimintestData$metadata[15:23])))


## ----richness_diet_run, fig.cap="Microbial richness as function of diet over time", fig.height=3, fig.width=9.7----

alimintestData$metadata$richness       = rep(NA,76)
alimintestData$metadata$richness = round(rarefy(t(alimintestData$otu), 1000))
richness_diet_run = ggplot(alimintestData$metadata, aes(x=Diet_run, y=richness)) + geom_boxplot(aes(fill=Diet_run)) + facet_grid(.~Time.point) + scale_fill_manual(values=c("orange","blue")) + theme_bw()
richness_diet_run


## ----rarefaction_analysis, fig.cap="Rarefaction analysis", fig.height=6, fig.width=9.7----

res  = round(rarefy(t(alimintestData$otu), seq(10,1000,10)))
colnames(res) = seq(10,1000,10)
res  = cbind(res,alimintestData$metadata[,1:5])
res2 = melt(res, id.vars=c(colnames(alimintestData$metadata[,1:5])))

res2$variable   = as.numeric(as.character(res2$variable))
res2$Time.point = paste('Time point',res2$Time.point, sep="")

ggplot(res2, aes(x=variable, y=value, group=Sample_ID)) +
 geom_line(aes(color=Diet_run)) + scale_color_manual("",values=c("orange","blue")) +
 facet_wrap(~ Time.point, ncol=2) + theme_bw() + xlab("Number of sequences sampled") + ylab("nb of OTUs")

 

## ----microbiota_qpcr, fig.cap="Gut microbiota composition assessed by qPCR as function of diet over time", fig.height=9.7, fig.width=9.7----

a = alimintestData$metadata[c(3,5,15:21)]
b = alimintestData$basalqpcr[c(3,5,6:12)]

colnames(a) = colnames(b)

df = rbind(a,b)


qPCR_names=c("log_Cleptum_group","log_Bcoccoides_group","log_Bacteroides_Prevotella","log_Bifidobacteria","log_Ecoli","log_Lactobacillus")

for(i in qPCR_names) {

  df[,i] = df[,i] - df$log_all_Bacteria
}

colnames(df)[4:9] = c("Cleptum_group","Bcoccoides_group","Bacteroides_Prevotella","Bifidobacteria","Ecoli","Lactobacillus")

df = melt(df, id.vars=c("Diet_run","Time.point"))

qPCR_raw_boxplot=ggplot(df, id.vars=c("Time.point","Diet_run")) + geom_boxplot(aes(y=value,x=Diet_run, fill=Diet_run)) + facet_grid(variable~Time.point, scales="free") + scale_fill_manual(values=c("orange","blue")) + theme_bw()

qPCR_raw_boxplot


## ----microbiota_qpcr_richness, fig.cap="Gut microbiota richness and composition assessed by 454 and qPCR as function of diet over time", fig.height=12.7, fig.width=9.7, fig.show="hide"----

grid.arrange(richness_diet_run + theme(legend.position = "none") + xlab(""), qPCR_raw_boxplot + xlab("") + theme(legend.position = "none") + ylab("qPCR assays"), heights=c(2/8, 6/8))



## ----microbiota_tax_bca, include=TRUE, echo=TRUE-------------------------

otu_norm = prop.table(as.matrix(alimintestData$otu),2)

microbiota_tax = list(
	phylum_tax = apply(otu_norm, 2, tapply, alimintestData$tax$phylum, sum ),
	class_tax  = apply(otu_norm, 2, tapply, alimintestData$tax$class, sum ),
	order_tax  = apply(otu_norm, 2, tapply, alimintestData$tax$order, sum ),
	family_tax = apply(otu_norm, 2, tapply, alimintestData$tax$family, sum ),
	genus_tax  = apply(otu_norm, 2, tapply, alimintestData$tax$genus, sum ),
	otu_norm   = noise.removal(otu_norm, percent=1)
)

metadata76            = alimintestData$metadata
metadata76$Time.point = as.factor(metadata76$Time.point)

microbiota_tax_variations=alimintest_metadata_pca=sapply(microbiota_tax, function(x) {

x                  = as.data.frame(t(log10(x + 10^-5)))
x.pca              = dudi.pca(x,    scannf=FALSE, nf=2)
x.bca.ind          = bca(x.pca,     scannf=FALSE, nf=2, fac=as.factor(metadata76$Subject_Id))
x.wca.ind          = bca(x.pca,     scannf=FALSE, nf=2, fac=as.factor(metadata76$Subject_Id))
x.bca.diet         = bca(x.pca,     scannf=FALSE, nf=2, fac=as.factor(metadata76$Diet_run:metadata76$Time.point))
#x.wca.ind.bca.diet = bca(x.wca.ind, scannf=FALSE, nf=2, fac=as.factor(metadata76$Diet_run:metadata76$Time.point))
x.wca.ind.bca.diet = bca(wca(dudi.pca(x, scannf=F, nf=2), scannf=F, nf=2, fac=as.factor(metadata76$Subject_Id)), scannf=F, nf=2, fac=as.factor(metadata76$Diet_run:metadata76$Time.point) )

x.bca.ind.res          = randtest(x.bca.ind) 
x.bca.diet.res         = randtest(x.bca.diet) 
x.wca.ind.bca.diet.res = randtest(x.wca.ind.bca.diet)

res = data.frame(ind.obs        = x.bca.ind.res$obs*100,
                 ind.sim        = quantile(x.bca.ind.res$sim, 0.95)[[1]]*100,
		         diet.obs       = x.bca.diet.res$obs*100,
		         diet.sim       = quantile(x.bca.diet.res$sim,0.95)[[1]]*100,
				 diet.w.ind.obs = x.wca.ind.bca.diet.res$obs*100,
				 diet.w.ind.sim = quantile(x.wca.ind.bca.diet.res$sim, 0.95)[[1]]*100
				 )
return(res)
}
)


## ----between_class_analysis_tax3, fig.cap="Between class analysis using diet and individual as factors", fig.height=5.5, fig.width=7----

m           = data.frame(matrix(unlist(microbiota_tax_variations), nrow = nrow(microbiota_tax_variations), byrow=F))
colnames(m) = colnames(microbiota_tax_variations)
rownames(m) = rownames(microbiota_tax_variations)

microbiota_tax_variations = m

m          = melt(data.frame(microbiota_tax_variations, id=rownames(microbiota_tax_variations)), id="id", value.name = "Inertia")
m$DataType = rep(c("observed","MC simulated\n(95th percentile)"),18)
m$Test     = rep(rep(c("Between Subject","Between Diet","Between Diet\nwithin Subject"),each=2),6)


ggplot(m, aes(x=variable, y=Inertia, shape=Test, col=DataType)) + geom_point(size=3, alpha=0.7) + xlab("") + ylab("Inertia attributed (%)") + scale_x_discrete(label=c("phylum","class","order","family","genus","OTU")) + theme_bw()


## ----coinertia_analysis--------------------------------------------------
genus.jsd = melt(as.matrix(BiotypeR::dist.JSD(microbiota_tax$genus_tax)))

microbial.change=NULL

metadata76$Subject_Id = as.factor(as.character(metadata76$Subject_Id))


for (s in levels(metadata76$Subject_Id)) {

		t2 = metadata76$Sample_ID[which(metadata76$Subject_Id==s & metadata76$Time.point=="2")]
		t3 = metadata76$Sample_ID[which(metadata76$Subject_Id==s & metadata76$Time.point=="4")]
		t4 = metadata76$Sample_ID[which(metadata76$Subject_Id==s & metadata76$Time.point=="3")]
		t5 = metadata76$Sample_ID[which(metadata76$Subject_Id==s & metadata76$Time.point=="5")]
		
		
		jsd=c(
		genus.jsd[which(genus.jsd[,1]==t2 & genus.jsd[,2]==t3),3],
		genus.jsd[which(genus.jsd[,1]==t3 & genus.jsd[,2]==t4),3],
		genus.jsd[which(genus.jsd[,1]==t4 & genus.jsd[,2]==t5),3]
		)
		
		richness.before=c(metadata76$richness[metadata76$Sample_ID==t2],
		metadata76$richness[metadata76$Sample_ID==t3],
		metadata76$richness[metadata76$Sample_ID==t4])
		
		richness.after=c(
		metadata76$richness[metadata76$Sample_ID==t3],
		metadata76$richness[metadata76$Sample_ID==t4],
		metadata76$richness[metadata76$Sample_ID==t5]
		)
		
		comet.after=c(
		metadata76$Comet_assay[metadata76$Sample_ID==t3],
		metadata76$Comet_assay[metadata76$Sample_ID==t4],
		metadata76$Comet_assay[metadata76$Sample_ID==t5]
		)
		
		
		diet=metadata76$Diet_run[which(metadata76$Subject_Id==s)][1]
		
		tmp=data.frame(Subject_ID=rep(s,3),change=c("2_3","3_4","4_5"),diet=rep(diet,3),jsd=jsd, richness.before, richness.after, comet.after)
		
		microbial.change=rbind(microbial.change,tmp)
		
}

METHOD="spearman"

cor.test(log10(microbial.change$jsd), microbial.change$richness.before, method=METHOD)
cor.test(log10(microbial.change$jsd), microbial.change$comet.after, method=METHOD)


## ----richness_vs_jsd, fig.cap="Microbial change vs OTU richness", fig.height=4, fig.width=14----
group         = as.factor(microbial.change$change:microbial.change$diet)
levels(group) = c("after 10g fiber/day","after 40g fiber/day", "usual diet (washout)", 
                  "usual diet (washout)", "after 40g fiber/day","after 10g fiber/day")

richness_vs_jsd = ggplot(data.frame(microbial.change,group), aes(x=log10(jsd), y = richness.before)) +
                  geom_point(aes(colour = diet), cex = 4) + stat_smooth(method = "lm") + 
                  scale_colour_manual(values = c("orange","blue")) +
                  facet_grid(.~group) + 
                  xlab("microbiota change\n(estimated by log10 JSD distance)") +
                  ylab("Microbial richness") +
				  theme_bw()

richness_vs_jsd


## ----comet_vs_jsd, fig.cap="Microbial change vs comet assays", fig.height=4, fig.width=14, fig.show="hide"----
ggplot(data.frame(microbial.change,group), aes(x=log10(jsd), y=comet.after)) +
  geom_point(aes(colour=diet), cex=4) + stat_smooth(method = "lm") +
  xlab("microbiota change\n(estimated by log10 JSD distance)")  +
  ylab("comet assays") + theme_bw() + facet_grid(.~group)

## ----cor_test_change, include=FALSE--------------------------------------

x1 = cor.test(log10(microbial.change$jsd[which(group=="after 10g fiber/day")]),  
microbial.change$richness.before[which(group=="after 10g fiber/day")], method=METHOD)
y1 = cor.test(log10(microbial.change$jsd[which(group=="after 40g fiber/day")]),  
microbial.change$richness.before[which(group=="after 40g fiber/day")], method=METHOD)
z1 = cor.test(log10(microbial.change$jsd[which(group=="usual diet (washout)")]),  
microbial.change$richness.before[which(group=="usual diet (washout)")], method=METHOD)

x2 = cor.test(log10(microbial.change$jsd[which(group=="after 10g fiber/day")]),  
microbial.change$comet.after[which(group=="after 10g fiber/day")], method=METHOD)
y2 = cor.test(log10(microbial.change$jsd[which(group=="after 40g fiber/day")]),  
microbial.change$comet.after[which(group=="after 40g fiber/day")], method=METHOD)
z2 = cor.test(log10(microbial.change$jsd[which(group=="usual diet (washout)")]),  
microbial.change$comet.after[which(group=="usual diet (washout)")], method=METHOD)

x1 = paste(round(x1$estimate[[1]],2), " (",round(x1$p.value,3),")", sep="")
y1 = paste(round(y1$estimate[[1]],2), " (",round(y1$p.value,3),")", sep="")
z1 = paste(round(z1$estimate[[1]],2), " (",round(z1$p.value,3),")", sep="")

x2 = paste(round(x2$estimate[[1]],2), " (",round(x2$p.value,3),")", sep="")
y2 = paste(round(y2$estimate[[1]],2), " (",round(y2$p.value,3),")", sep="")
z2 = paste(round(z2$estimate[[1]],2), " (",round(z2$p.value,3),")", sep="")

res = t(matrix(c(x1,x2,y1,y2,z1,z2), ncol=3))

colnames(res) = c("Richness before diet","comet assays after diet")
rownames(res) = c("10g","40g","usual")

res



## ----cor_test_change_kable, results="asis", comment=NA, echo=FALSE-------

kable(format(res))


## ----comet_richness_PC1--------------------------------------------------
genus.scfa.coi = coinertia(
		         dudi.pca(na.omit(metadata76[7:14]), scannf=F, nf=7), 
		         dudi.pca(as.data.frame(t(log10(noise.removal(microbiota_tax$genus_tax, percent=1)[,-attr(na.omit(metadata76[7:14]),"na.action")] + 10^-5))), scannf=F, nf=7),
		         scannf=F, nf=7
		)

test=randtest(genus.scfa.coi,9999)


dd=data.frame(genus.scfa.coi$lX, metadata76[-c(50,70),], as.data.frame(t(log10(noise.removal(microbiota_tax$genus_tax, percent=1)[,-attr(na.omit(metadata76[7:14]),"na.action")] + 10^-5))))


## ----richness_vs_SCFA_genus, fig.cap="association between OTU richness, SCFA and bacterial genera", fig.height=4, fig.width=13----

p1 = ggplot(dd, aes(x = AxcX1, y = richness)) +
     geom_point(aes(size = richness, colour = caproate+valerate), alpha = 0.8) +
     scale_colour_gradient("Caproate and\nValerate",low = "blue", high = "red", space = "Lab") +
     xlab("PC1") + ylab("Microbial richness") #+ xlim( -4,3) 
	 
p2 = ggplot(dd, aes(x = AxcX1, y = richness)) +
     geom_point(aes(size = richness, colour = Prevotella-Bacteroides), alpha = 0.8) +
     scale_colour_gradient("log10 of\nPrevotella/\nBacteroides\nratio",low = "blue", high = "red", space = "Lab") +
     xlab("PC1") + ylab("Microbial richness") #+ xlim( -4,3) 
	 
grid.arrange(p1+theme_bw(),p2+theme_bw(), ncol=2)


