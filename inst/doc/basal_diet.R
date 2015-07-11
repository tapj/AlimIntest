## ----load_library_data---------------------------------------------------
#analysis
library(ade4)
library(ggplot2)
library(reshape2)
library(gridExtra)
data(alimintestData)
data(basal_diet)

## ------------------------------------------------------------------------

metadata76            = alimintestData$metadata
metadata76$Time.point = as.factor(metadata76$Time.point)
basal_diet_tmp        = basal_diet[which(basal_diet$visite=="before_study"),]

for(i in 4:dim(basal_diet_tmp)[2]){

	basal_diet_tmp[,i] = as.factor(as.character(basal_diet_tmp[,i]))

}

diet.acm = dudi.acm(basal_diet_tmp[,-c(1:3)], nf=3, scannf=FALSE)
diet.bca = bca(diet.acm, fac=as.factor(basal_diet_tmp[,1]), nf=3, scannf=F)
#randtest(diet.bca)

basal_richness     = metadata76[which(metadata76$Time.point=="2"),c("Subject_Id","richness")]
diet.richness.bca  = merge(data.frame(Subject_Id=row.names(diet.bca$li), diet.bca$li), basal_richness, by="Subject_Id")
#diet.bca$co[order(diet.bca$co[,3]+diet.bca$co[,2]),]
diet.richness.bca2 = merge(data.frame(diet.bca$ls, Subject_Id=basal_diet_tmp[,1]),diet.richness.bca, by= "Subject_Id")


## ----basal_diet_richness, fig.cap="Associations between basal diet and microbiota richness", fig.height=7, fig.width=14----

p1_diet = ggplot(diet.richness.bca2) + 
geom_point(aes(x=CS2,y=CS3)) + 
geom_segment(aes(x=CS2,y=CS3, xend=Axis2, yend=Axis3, col=Subject_Id)) + 
geom_text(aes(x=Axis2, y=Axis3, label=Subject_Id)) + xlab("PC2") +ylab("PC3") + guides(col=FALSE)

p2_diet = ggplot(diet.richness.bca, aes(x=Axis2+Axis3, y=richness, col=richness)) + 
geom_point(size=4) + xlab("PC2 + PC3") + ylab("Microbiota richness") 

grid.arrange(p1_diet,p2_diet, ncol=2)


