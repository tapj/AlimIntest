###############################################################
#                                                             #
#           Command Alimintest paper for R                    #
#                                                             #
###############################################################

library(ade4)



## I - Data preparation
#read.csv2("table-global-alimintest-may2010.csv")-> alimintest
read.csv2("table-global-alimintest-june2010.csv")-> alimintest
read.csv2("methods.csv")-> methods

row.names(alimintest)<-alimintest$SampleID
alimintest$subject<-factor(alimintest$subject)
alimintest$diet<-factor(alimintest$diet)
alimintest$SampleID<-as.character(alimintest$SampleID)
alimintest$time<-factor(alimintest$time)


alimintest[,methods$type=="qPCR"] -> alimintest.qpcr
alimintest[,methods$type=="SCFA_fecal"] -> alimintest.SCFA_fecal
alimintest[,methods$type=="SCFA_water"] -> alimintest.SCFA_water
alimintest[,methods$type=="454_classifier"] -> alimintest.genus
alimintest[,methods$dominant=="yes"] -> alimintest.dominant
alimintest[,methods$type=="454_core"] -> alimintest.core
alimintest[,methods$type=="model_euk"] -> alimintest.model



## II - Co-inertia analysis : compute inertia shared between two methods

TableS2 <- matrix(nr=7, nc=7)

# qPCR vs SCFA_fecal

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[2,1]

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,2]

# qPCR vs SCFA_water

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[3,1]

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,3]

# qPCR vs alimintest.genus

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[4,1]

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,4]

# qPCR vs alimintest.dominant

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[5,1]

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,5]

# qPCR vs alimintest.core

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[6,1]

randtest(coinertia(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$X454=="1"  ,2:6] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$qPCR=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,6]

# SCFA_fecal vs SCFA_water

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$SCFA_water=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$SCFA_fecal=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue  -> TableS2[3,2]

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$SCFA_water=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$SCFA_fecal=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs  -> TableS2[2,3]

# SCFA_fecal vs alimintest.genus

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[4,2]

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[2,4]

# SCFA_fecal vs alimintest.dominant

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[5,2]

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[2,5]

# SCFA_fecal vs alimintest.core

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[6,2]

randtest(coinertia(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[2,6]

# SCFA_water vs alimintest.genus

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[4,3]

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[3,4]

# SCFA_water vs alimintest.dominant

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[5,3]

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[3,5]

# SCFA_water vs alimintest.core

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[6,3]

randtest(coinertia(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[3,6]

# alimintest.genus vs alimintest.dominant

randtest(coinertia(dudi.pca(alimintest.genus[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[5,4]


randtest(coinertia(dudi.pca(alimintest.genus[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[4,5]

# alimintest.genus vs alimintest.core

randtest(coinertia(dudi.pca(alimintest.genus[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[6,4]

randtest(coinertia(dudi.pca(alimintest.genus[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[4,6]

# alimintest.dominant vs alimintest.core

randtest(coinertia(dudi.pca(alimintest.dominant[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[6,5]

randtest(coinertia(dudi.pca(alimintest.dominant[alimintest$X454=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$X454=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[5,6]

#alimintest.model vs qPCR

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$qPCR=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.qpcr[alimintest$model=="1" &  alimintest$qPCR=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,1]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$qPCR=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.qpcr[alimintest$model=="1" &  alimintest$qPCR=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[1,7]



#alimintest.model vs SCFA raw fecal

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$SCFA_fecal=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$model=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,2]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$SCFA_fecal=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$model=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[2,7]



#alimintest.model vs SCFA fecal water

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$SCFA_water=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$model=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,3]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$SCFA_water=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$model=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[3,7]



#alimintest.model vs total genera

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,4]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.genus[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[4,7]



#alimintest.model vs dominant

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,5]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.dominant[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[5,7]


#alimintest.model vs core species

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$pvalue -> TableS2[7,6]

randtest(coinertia(dudi.pca(alimintest.model[alimintest$model=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.core[alimintest$model=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))$obs -> TableS2[6,7]



# TableS2 creation
as.data.frame(TableS2)-> TableS2
row.names(TableS2)<-c("qPCR", "SCFA raw fecal", "SCFA fecal water", "Total genera", "Dominant genera", "Core species", "Model cells")
names(TableS2)<-c("qPCR", "SCFA raw fecal", "SCFA fecal water", "Total genera", "Dominant genera", "Core species", "Model cells")

write.csv2(TableS2, file="TableS2.csv")



## III - between and within class analysis

TableS1 <- matrix(nr=7, nc=6)


# qPCR
randtest(between(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$qPCR=="1"] , scannf=F, nf=3))$obs -> TableS1[1,1]

randtest(between(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$qPCR=="1"] , scannf=F, nf=3))$pvalue -> TableS1[1,2]

randtest(between(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$qPCR=="1"]:alimintest$diet[alimintest$qPCR=="1"]) , scannf=F, nf=3))$obs -> TableS1[1,3]


randtest(between(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$qPCR=="1"]:alimintest$diet[alimintest$qPCR=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[1,4]



randtest(between(within(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$qPCR=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$qPCR=="1"]:alimintest$diet[alimintest$qPCR=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[1,5]


randtest(between(within(dudi.pca(alimintest.qpcr[alimintest$qPCR=="1" ,1]-alimintest.qpcr[alimintest$qPCR=="1"  ,2:6] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$qPCR=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$qPCR=="1"]:alimintest$diet[alimintest$qPCR=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[1,6]



# SCFA raw fecal

randtest(between(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_fecal=="1"] , scannf=F, nf=3))$obs -> TableS1[2,1]


randtest(between(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_fecal=="1"] , scannf=F, nf=3))$pvalue -> TableS1[2,2]



randtest(between(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$SCFA_fecal=="1"]:alimintest$diet[alimintest$SCFA_fecal=="1"]) , scannf=F, nf=3))$obs -> TableS1[2,3]


randtest(between(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$SCFA_fecal=="1"]:alimintest$diet[alimintest$SCFA_fecal=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[2,4]




randtest(between(within(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_fecal=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$SCFA_fecal=="1"]:alimintest$diet[alimintest$SCFA_fecal=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[2,5]


randtest(between(within(dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_fecal=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$SCFA_fecal=="1"]:alimintest$diet[alimintest$SCFA_fecal=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[2,6]





# SCFA fecal water


randtest(between(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$SCFA_water=="1"])), scannf=F, nf=3))$obs -> TableS1[3,1]


randtest(between(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$SCFA_water=="1"])), scannf=F, nf=3))$pvalue -> TableS1[3,2]


randtest(between(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$SCFA_water=="1"]:alimintest$diet[alimintest$SCFA_water=="1"]) , scannf=F, nf=3))$obs -> TableS1[3,3]

randtest(between(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$SCFA_water=="1"]:alimintest$diet[alimintest$SCFA_water=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[3,4]

randtest(between(within(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_water=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$SCFA_water=="1"]:alimintest$diet[alimintest$SCFA_water=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[3,5]

randtest(between(within(dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$SCFA_water=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$SCFA_water=="1"]:alimintest$diet[alimintest$SCFA_water=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[3,6]

# total genera

randtest(between(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$obs -> TableS1[4,1]



randtest(between(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$pvalue -> TableS1[4,2]



randtest(between(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$obs -> TableS1[4,3]


randtest(between(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[4,4]




randtest(between(within(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[4,5]

randtest(between(within(dudi.pca(alimintest.genus[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[4,6]





# dominant genera

randtest(between(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$obs -> TableS1[5,1]

randtest(between(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$pvalue -> TableS1[5,2]


randtest(between(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$obs -> TableS1[5,3]


randtest(between(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[5,4]


randtest(between(within(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[5,5]

randtest(between(within(dudi.pca(alimintest.dominant[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[5,6]





# dominant core

randtest(between(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$obs -> TableS1[6,1]



randtest(between(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$X454=="1"])), scannf=F, nf=3))$pvalue -> TableS1[6,2]



randtest(between(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$obs -> TableS1[6,3]



randtest(between(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[6,4]



randtest(between(within(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[6,5]


randtest(between(within(dudi.pca(alimintest.core[alimintest$X454=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$X454=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$X454=="1"]:alimintest$diet[alimintest$X454=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[6,6]




# Model eukaryotic cells

randtest(between(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$model=="1"])), scannf=F, nf=3))$obs -> TableS1[7,1]



randtest(between(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=factor(as.vector(alimintest$subject[alimintest$model=="1"])), scannf=F, nf=3))$pvalue -> TableS1[7,2]



randtest(between(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$model=="1"]:alimintest$diet[alimintest$model=="1"]) , scannf=F, nf=3))$obs -> TableS1[7,3]



randtest(between(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=factor(alimintest$time[alimintest$model=="1"]:alimintest$diet[alimintest$model=="1"]) , scannf=F, nf=3))$pvalue -> TableS1[7,4]



randtest(between(within(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$model=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$model=="1"]:alimintest$diet[alimintest$model=="1"]) ,  scannf=F, nf=3))$obs -> TableS1[7,5]


randtest(between(within(dudi.pca(alimintest.model[alimintest$model=="1" ,] , scannf=F, nf=3) , fac=alimintest$subject[alimintest$model=="1"] , scannf=F, nf=3), fac=factor(alimintest$time[alimintest$model=="1"]:alimintest$diet[alimintest$model=="1"]) ,  scannf=F, nf=3))$pvalue -> TableS1[7,6]




as.data.frame(TableS1)-> TableS1

row.names(TableS1)<-c("qPCR", "SCFA raw fecal", "SCFA fecal water", "Total genera", "Dominant genera", "Core species", "Model euk cells")

names(TableS1)<-c("Between subject","pvalue", "Between diet", "pvalue", "Between diet within subject", "pvalue")


write.csv2(TableS1, "TableS1.csv")



### Figure 2 qPCR between class analysis
### TO PUT HERE


### Figure 3 PTA

names(alimintest.dominant)<- c("Bacteroides","Bifidobacterium", "Coprobacillus", "Blautia", "Dorea","Roseburia", "Parabacteroides", "Prevotella", "Alistipes", "Faecalibacterium", "Oscillibacter", "Ruminococcus",  "Subdoligranulum","Dialister" )





withinpca(data.frame(alimintest.dominant[20:95,]), factor(as.character(alimintest$time[20:95])), scal="partial", scannf=F)-> diet.time.pca

#withinpca(data.frame(alimintest.dominant[20:95,]), factor(factor(as.character(alimintest$time[20:95])):factor(as.character(alimintest$diet_desc[20:95]))), scal="partial", scannf=F)-> diet.time.pca #combine time and diet run

diet.time.kta <- ktab.within(diet.time.pca, colname=rep(1:19,4)  )



taille=0.75
pdf("PTA-microbiote-in-time_may2010.pdf", h=10, w=12, v="1.4")
plot(pta(t(diet.time.kta), scannf=F, nf=3)$Tli[,1], pta(t(diet.time.kta), scannf=F, nf=3)$Tli[,2], type='n', xlab="PC1 21.38 %", ylab="PC2 16.16 %", xlim=c(-4.5,4), ylim=c(-3,4), font.lab=2)

s.class(pta(t(diet.time.kta), scannf=F, nf=3)$Tli, fac=factor(as.character(alimintest$diet[20:95])), cell=0, cstar=0, cpoi=2,  clab=0, col=c("blue", "orange"), grid = F, add.plot=T)


s.match(pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F, nf=3)$TL[,1]==1,], pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F, nf=3)$TL[,1]==2,], edge=F, clab=taille, label=paste(c(1,3:20),".D1"), xax=1, yax=2, add.plot=T , cpoi=0)


s.match(pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F, nf=3)$TL[,1]==2,], pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F, nf=3)$TL[,1]==3,], add.plot=T, edge=F, clab=taille,label=paste(c(1,3:20),".W"), xax=1, yax=2, cpoi=0)

s.match(pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F)$TL[,1]==3,], pta(t(diet.time.kta), scannf=F, nf=3)$Tli[pta(t(diet.time.kta), scannf=F, nf=3)$TL[,1]==4,], add.plot=T, edge=T, clab=taille, label=paste(c(1,3:20),".D2"), xax=1,yax=2, cpoi=0)


text(pta(t(diet.time.kta), scannf=F, nf=3)$co[,1]*2.5, pta(t(diet.time.kta), scannf=F, nf=3)$co[,2]*2.5, labels=row.names(pta(t(diet.time.kta), scannf=F, nf=3)$co), font=4, pos=2)

legend(-4,4, c("D1: 1st diet period", "W: Wash-out period", "D2: 2nd diet period", "40-10 sequence", "10-40 sequence"), bty='o', pch=16, pt.cex=c(0,0,0,2,2), col=c( "black","black", "black", "blue", "orange") )



dev.off()


### Figure 4 CorCircle !!

pdf("Figure4_SCFA_microbiota_may2010.pdf", h=12, w=12, v="1.4")
par(mfrow=c(2,2))
s.corcircle(coinertia(dudi.pca(-(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,2:6]) , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1, grid=F)


s.corcircle(coinertia(dudi.pca(-(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1"  ,2:6]) , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$qPCR=="1" &  alimintest$SCFA_fecal=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1, add.plot=T)


s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1, grid=F)


s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_fecal[alimintest$SCFA_fecal=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1, add.plot=T)

#dev.off()



#pdf("Figure4_SCFA-w_microbiota_may2010.pdf", h=6, w=12, v="1.4")
#par(mfrow=c(1,2))
s.corcircle(coinertia(dudi.pca(-(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,2:6]) , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1, grid=F)


s.corcircle(coinertia(dudi.pca(-(alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,1]-alimintest.qpcr[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1"  ,2:6]) , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$qPCR=="1" &  alimintest$SCFA_water=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1, add.plot=T)


s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1, grid=F)


s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$SCFA_water=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.SCFA_water[alimintest$SCFA_water=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1, add.plot=T)

dev.off()







### Figure S2 Enterostates
##TO PUT HERE

### Figure 6 Test comet

test<-matrix(nr=210,nc=4)

for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" | methods$type=="qPCR")) {cor(alimintest$comet[alimintest$SCFA_water=="1"], alimintest[alimintest$SCFA_water=="1", i], method="spearman") -> test[i,1]; cor.test(alimintest$comet[alimintest$SCFA_water=="1"], alimintest[alimintest$SCFA_water=="1", i], method="spearman")$p.value-> test[i,2] ;  wilcox.test(alimintest[alimintest$SCFA_water=="1",i]~alimintest$comet[alimintest$SCFA_water=="1"]>50)$p.value -> test[i,3]; names(alimintest[i])-> test[i,4]  }


pdf("Figure_Supp_comet.pdf", h=14, w=8)
par(mfrow=c(4,2))

for ( i in which(test[,3]<0.05 & test[,2]<0.05)) 
{stripchart(alimintest[alimintest$SCFA_water=="1",i]~alimintest$comet[alimintest$SCFA_water=="1"]>50, method="jitter", vert=T, pch=17, ylab=names(alimintest[i]), xlab="comet > 50%" , cex=1.5)}



lm(alimintest$comet[alimintest$SCFA_water=="1"]~alimintest$UUAA2DE021[alimintest$SCFA_water=="1"]+alimintest$UUAB3DA091[alimintest$SCFA_water=="1"]+alimintest$UUAB3DA091[alimintest$SCFA_water=="1"]+alimintest$UUAP4BG111[alimintest$SCFA_water=="1"]  +alimintest$UUAV1AD121[alimintest$SCFA_water=="1"]  +alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water=="1"]   +alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water=="1"]   +alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water=="1"])



plot( alimintest$comet[alimintest$SCFA_water=="1"], 443.10*alimintest$UUAP4BG111[alimintest$SCFA_water == "1"] -177.60*alimintest$UUAA2DE021[alimintest$SCFA_water == "1"] -13.63*alimintest$UUAB3DA091[alimintest$SCFA_water == "1"] -135.38*alimintest$UUAV1AD121[alimintest$SCFA_water == "1"]+7075.25*alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1"]-3430.70*alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1"]+7894.33*alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1"] +43.45, ylab="expected comet % (7  taxa based model)" , xlab="observed comet %", pch=17, xlim=c(0,100), ylim=c(0,100))


dev.off()


## Figure Publi + photos

pdf("Figure6_comet_microbiota_June2010.pdf", h=6, w=11)

#layout(matrix(c(1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,2,2,3,3,4,4,5,5,6,6,2,2,3,3,4,4,5,5,6,6), 6,10, byrow=T))

#layout(matrix(c(1,1,0,0,0,1,1,0,0,0,2,3,4,5,6), 3,5, byrow=T))


layout(matrix(c(1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,2,2,3,3,4,4,5,5,6,6,7,7,8,8,2,2,3,3,4,4,5,5,6,6,7,7,8,8), 6,14, byrow=T))

par(mar=c(5,5,1,2))

plot( alimintest$comet[alimintest$SCFA_water=="1"], 443.10*alimintest$UUAP4BG111[alimintest$SCFA_water == "1"] -177.60*alimintest$UUAA2DE021[alimintest$SCFA_water == "1"] -13.63*alimintest$UUAB3DA091[alimintest$SCFA_water == "1"] -135.38*alimintest$UUAV1AD121[alimintest$SCFA_water == "1"]+7075.25*alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1"]-3430.70*alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1"]+7894.33*alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1"] +43.45, ylab="expected comet % (7  taxa based model)" , xlab="observed comet %", pch=17, xlim=c(0,100), ylim=c(-5,100),cex.lab=1.5, cex.axis=1.5)

ID=18
points(alimintest$comet[alimintest$SCFA_water=="1" & alimintest$SampleID == ID], 443.10*alimintest$UUAP4BG111[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -177.60*alimintest$UUAA2DE021[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -13.63*alimintest$UUAB3DA091[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -135.38*alimintest$UUAV1AD121[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]+7075.25*alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]-3430.70*alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]+7894.33*alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] +43.45, pch=1, cex=3, lwd=2, col="blue")




ID=188

points(alimintest$comet[alimintest$SCFA_water=="1" & alimintest$SampleID == ID], 443.10*alimintest$UUAP4BG111[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -177.60*alimintest$UUAA2DE021[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -13.63*alimintest$UUAB3DA091[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] -135.38*alimintest$UUAV1AD121[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]+7075.25*alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]-3430.70*alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID]+7894.33*alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID] +43.45, pch=1, cex=3, lwd=2, col="red")

ID1=18
ID2=188

par(mar=c(5,5,2,2))

barplot(c(alimintest$UUAB3DA091[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100, alimintest$UUAB3DA091[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), ylab="% abundance", xlab="F. prausnitzii",cex.lab=1.5, cex.axis=1.5)

barplot(c(alimintest$UUAA2DE021[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100,alimintest$UUAA2DE021[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), xlab="B. longum",cex.lab=1.5, cex.axis=1.5)


barplot(c(alimintest$UUAV1AD121[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100,alimintest$UUAV1AD121[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), xlab="E. eligens",cex.lab=1.5, cex.axis=1.5)


barplot(c(alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100,alimintest$Erysipelotrichaceae_Turicibacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), xlab="Turicibacter",cex.lab=1.5, cex.axis=1.5)


barplot(c(alimintest$UUAP4BG111[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100,alimintest$UUAP4BG111[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), xlab="Ruminococcus sp.",cex.lab=1.5, cex.axis=1.5)


barplot(c(alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100 ,alimintest$Ruminococcaceae_Sporobacter[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100 ), col=c("blue", "red"), xlab="Sporobacter",cex.lab=1.5, cex.axis=1.5)


barplot(c(alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1" & alimintest$SampleID == ID1]*100,alimintest$Desulfovibrionaceae_Desulfovibrio[alimintest$SCFA_water == "1" & alimintest$SampleID == ID2]*100), col=c("blue", "red"), xlab="Desulfovibrio",cex.lab=1.5, cex.axis=1.5)

dev.off()




##############################################################
# Figure S3
# Microbiota and model cells gene rapporter relation
#
##############################################################



randtest(coinertia(dudi.pca(alimintest.dominant[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))

par(mfrow=c(1,2))

s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1)


s.corcircle(coinertia(dudi.pca(alimintest.dominant[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1)



cor(alimintest.dominant[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,], method="spearman" )



randtest(coinertia(dudi.pca(alimintest.core[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3))

par(mfrow=c(1,2))

s.corcircle(coinertia(dudi.pca(alimintest.core[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$c1)


s.corcircle(coinertia(dudi.pca(alimintest.core[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , scannf=F, nf=3) , dudi.pca(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,] , scannf=F, nf=3) , scannf=F, nf=3)$l1)



cor(alimintest.core[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,], method="spearman" )



cor(alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1"  ,] , alimintest.model[alimintest$model_euk=="1" &  alimintest$X454=="1" ,], method="spearman" )




test2<-matrix(nr=211,nc=13)


#AP1
for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$AP1[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,2]; cor.test(alimintest$AP1[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,3] ;  names(alimintest[i])-> test2[i,1]  }

#FIAF
for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$FIAF[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,4]; cor.test(alimintest$FIAF[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,5] }

#NUR77
for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$NUR77[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,6]; cor.test(alimintest$NUR77[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,7] }

 
#NFkB
for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$NF.kB[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,8]; cor.test(alimintest$NF.kB[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,9]   }

#TSLP
for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$TSLP[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,10]; cor.test(alimintest$TSLP[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,11]   }


#PPARg

for ( i in which(methods$type=="454_classifier" | methods$type=="454_core" )) {cor(alimintest$PPARg[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1"& alimintest$X454=="1" , i], method="spearman") -> test2[i,12]; cor.test(alimintest$PPARg[alimintest$model_euk=="1" & alimintest$X454=="1" ], alimintest[alimintest$model_euk=="1" & alimintest$X454=="1" , i], method="spearman")$p.value-> test2[i,13]   }



data.frame(na.omit(test2))-> test2


names(test2)<- c("taxa_ID","AP1","AP1p","FIAF","FIAFp", "NUR77", "NUR77p", "NFkB", "NFkBp", "TSLP", "TSLPp", "PPARg", "PPARgp")



write.csv(file="model.euk.test.csv", data.frame( row.names=test2[,1], AP1=as.vector(test2$AP1), FIAF=as.vector(test2$FIAF), NUR77=as.vector(test2$NUR77), NFkB=as.vector(test2$NFkB),  TSLP=as.vector(test2$TSLP), PPARg=as.vector(test2$PPARg),  

AP1p=as.vector(test2$AP1p), FIAFp=as.vector(test2$FIAFp), NUR77p=as.vector(test2$NUR77p), NFkBp=as.vector(test2$NFkBp),  TSLPp=as.vector(test2$TSLPp), PPARgp=as.vector(test2$PPARgp)  

))

read.csv("model.euk.test.csv")-> model.euk.test

row.names(model.euk.test) <- model.euk.test[,1]
model.euk.test[-1]-> model.euk.test


p<-0.01
c<-0.3

model.euk.test[model.euk.test$AP1p< p | model.euk.test$FIAFp< p | model.euk.test$NUR77p< p | model.euk.test$NFkBp< p | model.euk.test$TSLPp< p |  model.euk.test$PPARgp< p , ] ->model.euk.test


model.euk.test[model.euk.test$AP1< -c | model.euk.test$FIAF< -c | model.euk.test$NUR77< -c | model.euk.test$NFkB< -c | model.euk.test$TSLP< -c |  model.euk.test$PPARg< -c  | model.euk.test$AP1 > c | model.euk.test$FIAF > c | model.euk.test$NUR77 > c | model.euk.test$NFkB > c | model.euk.test$TSLP > c |  model.euk.test$PPARg > c , ] ->model.euk.test


heatmap(as.matrix(model.euk.test[1:6]), scale="none")

##############################################################
# Figure S4 and table S1
#  rarefaction and eco diversity table
#
##############################################################

load("classifier.RData")
row.names(classifier)=sort(row.names(alimintest[alimintest$X454=="1",methods$type=="454_classifier"]))

classifier.metadata=alimintest[row.names(classifier),3:5]

par(mfrow=c(2,2))

coldiet=c("blue","orange")

for ( tp in 2:5) {

plot(1:10,1:10, xlim=c(1,12081), ylim=c(1,53), type="n", xlab="Number of sequences sampled", ylab="Number of genera", cex.lab=1.5)
legend("topright", paste("time point",tp), bty="n", cex=1.5)
legend("bottomright", c("diet 10-40","diet 40-10"), col=c("orange","blue"), lwd=1, bty="n", cex=1.5)


for (i in 1:dim(classifier[classifier.metadata$time==tp,])[1]) {



lines(as.vector(rarefy(data.frame(classifier[classifier.metadata$time==tp,][i,]), 1:sum(classifier[classifier.metadata$time==tp,][i,]), MARGIN=2)), col=coldiet[classifier.metadata[classifier.metadata$time==tp,2][i]]) 



}

}





require(vegan)

sample<-apply(t(classifier),1,sum)
#sample[6]<-100000

resolution=1000

rS2<-matrix(nr=round(max(sample)/resolution, 0)+1, nc=76)
rSe2<-matrix(nr=round(max(sample)/resolution, 0)+1, nc=76)
sample.matrix<-matrix(nr=round(max(sample)/resolution, 0)+1, nc=76)


pdf(file="alimintest_rarefaction.pdf", h=8, w=12)
 
plot( 1:10, 1:10, xlim=c(1,1000) , ylim=c(1,max(specnumber(t(classifier)))), type="n", xlab="Number of sequences sampled", ylab="Number of genus")

for (i in 1:76) {
r<-NULL

rarefy(data.frame(t(classifier))[i,], c(1,(1:round(sample[i]/resolution,0))*resolution), se=T)-> r
r[1,]->rS2[1:length(r[1,]),i]
r[2,]->rSe2[1:length(r[2,]),i]

lines(c(1,(1:round(sample[i]/resolution,0))*resolution),na.omit(rS2[,i]))

}

dev.off()





data.frame(Sample_ID=row.names(classifier), sequence_number=apply(t(classifier),1,sum), detected_genus=specnumber(classifier) ,Simpson=round(diversity(classifier, "simpson"), 3), Shannon=round(diversity(classifier), 3), Chao1=estimateR(classifier[1:76,])[2,]) -> ecoalimintest



write.csv2(file="ecoalimintest1.csv", ecoalimintest)


# Figure density

names(alimintest.dominant)<- c("Bacteroides","Bifidobacterium", "Coprobacillus", "Blautia", "Dorea","Roseburia", "Parabacteroides", "Prevotella", "Alistipes", "Faecalibacterium", "Oscillibacter", "Ruminococcus",  "Subdoligranulum","Dialister" )

par( mfrow=c(5,3))

for ( genus in names(alimintest.dominant)) {

genus.density=density(log10(na.omit(alimintest.dominant[,genus])))
plot(genus.density, xlab=paste(genus,"log10 proportion"), main="", ylab="Density")


}


# Figure replicats
pdf("replicats_2.pdf", h=7, w=14)

par(mfrow=c(1,2))


read.csv2("classifier_prop_cogenics_genoscreen.csv")->classifier2

classifier2[ classifier2[,2] >= 0.01 |classifier2[,3] >= 0.01  | classifier2[,4] >= 0.01 | classifier2[,5] >= 0.01 |classifier2[,6] >= 0.01|classifier2[,7] >= 0.01 | classifier2[,8] >= 0.01 | classifier2[,9] >= 0.01   , ] ->  classifier2


plot( log10(c(classifier2[,2],classifier2[,4], classifier2[,6], classifier2[,8],classifier2[,10],classifier2[,12], classifier2[,14], classifier2[,16])) , log10(c(classifier2[,3],classifier2[,5], classifier2[,7], classifier2[,9], classifier2[,11],classifier2[,13], classifier2[,15], classifier2[,17])), xlab="454 GS FLX Ti Compagny A", ylab="454 GS FLX Ti Compagny B", cex.lab=1.5, type="n", xlim=c(-4,0), ylim=c(-4,0))

test<-matrix(nr=8, nc=2)
for(i in 1:8)
{points(log10(classifier2[,i*2]),log10(classifier2[,(i*2)+1]), col=1, pch=16)
cor(classifier2[,i*2],classifier2[,(i*2)+1],  method="spearman")-> test[i,1]
cor.test(classifier2[,i*2],classifier2[,(i*2)+1],  method="spearman")$p.value -> test[i,2]
}

abline(0,1, lty=2)
mtext("a.", 3, at=-4.5, adj=1, font=2, cex=2)

read.csv2("classifier_prop_sanger_454.csv")-> classifier
classifier[classifier[,2] >= 0.01 |  classifier[,4] >= 0.01 | classifier[,5] >= 0.01 |classifier[,6] >= 0.01, c(1,2,4,5,6)] ->  classifier
cor(classifier[,2:5], method="spearman")


plot(log10(classifier$OF.Genoscope.Sanger), log10(classifier$OF.Genoscreen.Run1a), type="n", xlab="Sanger method", ylab="454 GS FLX Ti method", cex.lab=1.5, xlim=c(-2,-0.5), ylim=c(-2,-0.5))



points(log10(classifier$OF.Genoscope.Sanger), log10(classifier$OF.Genoscreen.Run1b), col="green", pch=16, cex=1.5)
points(log10(classifier$OF.Genoscope.Sanger), log10(classifier$OF.Genoscreen.Run2), col="blue", pch=16,cex=1.5)
legend("topleft", c("1/8 lanes with 12 MID", "1/8 lanes with 24 MID" , "1/16 lanes with 12 MID"), col=c("red","green","blue"), pch=16, bty='n', pt.cex=1.5)
points(log10(classifier$OF.Genoscope.Sanger), log10(classifier$OF.Genoscreen.Run1a), col="red", pch=16, cex=1.5)

abline(0,1, lty=2)
mtext("b.", 3, at=-2.2, adj=1, font=2, cex=2)

dev.off()












