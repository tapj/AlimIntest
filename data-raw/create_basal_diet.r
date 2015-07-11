basal_diet = read.csv2("data-raw/diet_basal_alimintest.csv")


basal_diet$sujet = as.character(basal_diet$sujet)
basal_diet$jour = as.character(basal_diet$jour)


for(i in 4:dim(basal_diet)[2]){

basal_diet[,i] = as.factor(basal_diet[,i])


}


save(basal_diet,file="data/basal_diet.RData")

