basal_diet = read.csv2("data-raw/diet_basal_alimintest.csv")


basal_diet$sujet = as.character(basal_diet$sujet)
basal_diet$jour = as.character(basal_diet$jour)


for(i in 4:dim(basal_diet)[2]){

basal_diet[,i] = as.factor(basal_diet[,i])


}


save(basal_diet,file="data/basal_diet.RData")



#analysis
library(ade4)
library(ggplot2)

basal_diet_tmp = basal_diet[which(basal_diet$visite=="before_study"),]


for(i in 4:dim(basal_diet_tmp)[2]){

basal_diet_tmp[,i] = as.factor(as.character(basal_diet_tmp[,i]))


}



diet.acm = dudi.acm(basal_diet_tmp[,-c(1:3)], nf=3, scannf=FALSE)
diet.bca = bca(diet.acm, fac=as.factor(basal_diet_tmp[,1]), nf=3, scannf=F)
randtest(diet.bca)


basal_richness = metadata76[which(metadata76$Time.point=="2"),c("Subject_Id","richness")]

diet.richness.bca = merge(data.frame(Subject_Id=row.names(diet.bca$li), diet.bca$li), basal_richness, by="Subject_Id")

ggplot(diet.richness.bca, aes(x=Axis1, y=Axis2, col=richness, size=richness)) + geom_point()

#test = merge(basal_richness,basal_diet[which(basal_diet$visite=="before_study" & basal_diet$jour=="-1"),], by.x="Subject_Id", by.y="sujet")

test_combine = apply(basal_diet[which(basal_diet$visite=="before_study"),-c(1:3)] , 2, function(x){ tapply(as.numeric(x), as.factor(basal_diet[which(basal_diet$visite=="before_study"),1]) ,sum) } )

test = merge(basal_richness, data.frame(Subject_Id=row.names(test_combine),test_combine)  )


#test = merge(basal_richness,basal_diet[which(basal_diet$visite=="before_study"),], by.x="Subject_Id", by.y="sujet")
#test = merge(basal_richness,basal_diet[which(basal_diet$visite=="after_study"),], by.x="Subject_Id", by.y="sujet")
#test = merge(basal_richness,basal_diet, by.x="Subject_Id", by.y="sujet")


for(i in 5:dim(test)[2]){

test[,i] = as.numeric(as.character(test[,i]))


}


test = test[, c(2,4+which(apply(test[-c(1:4)],2,sum)!=0))]


cor(test$richness, test[,-1], method="spearman")

apply(test[,-1], 2, function(x) cor.test(test$richness, x, method="spearman")$p.value)


diet_test_result = data.frame(rho=t(cor(test$richness, test[,-1], method="spearman")), p=apply(test[,-1], 2, function(x) wilcox.test(test$richness~x)$p.value))
diet_test_result = data.frame(rho=t(cor(test$richness, test[,-1], method="spearman")), p=apply(test[,-1], 2, function(x) cor.test(test$richness,x, method="spearman")$p.value))

diet_test_result[which(diet_test_result$p<0.05),]


ggplot(melt(test, id.vars=c("richness")), aes(fill=as.character(value), y=richness, x=variable)) + geom_boxplot() + coord_flip()


ggplot(data.frame(Richness = test$richness, dudi.pca(test[,-1], scannf=F, nf=10)$li)) + 
coord_tern() + geom_point(aes(
                  x=(Axis1-min(Axis1))/(max(Axis1)-min(Axis1)), 
                  y=(Axis2-min(Axis2))/(max(Axis2)-min(Axis2)), 
                  z=(Axis3-min(Axis3))/(max(Axis3)-min(Axis3)),  
                  size=test.richness)) + theme_rgbw() + 
xlab("PC1") + ylab("PC2") + zlab("PC3")


