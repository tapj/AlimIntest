# Generate metatrans_alimintest from raw data


# tax_metatrans = read_excel("~/storage/Dropbox/GitHub/AlimIntest/data-raw/metatranscriptomic_associated_data.xlsx", sheet =1)
# tax_16S = read_excel("~/storage/Dropbox/GitHub/AlimIntest/data-raw/metatranscriptomic_associated_data.xlsx", sheet = 2)
# metadata = read_excel("~/storage/Dropbox/GitHub/AlimIntest/data-raw/metatranscriptomic_associated_data.xlsx", sheet = 3)
 
tax_metatrans = read.csv2("data-raw/metatrans_taxonomy.csv")
tax_16S       = read.csv2("data-raw/metatrans_16S.csv")
metadata      = read.csv2("data-raw/metatrans_metadata.csv")
SubCatKegg    = prop.table(as.matrix(read.csv2("data-raw/SubCatKegg.csv", row=10)[2:9]),2)


tax_metatrans_summary = t(prop.table(apply(tax_metatrans[-which(tax_metatrans$manual_tax=="Discard"), 2:9], 2, tapply, as.character(tax_metatrans[-which(tax_metatrans$manual_tax=="Discard"),10]), sum),2))
tax_metatrans_summary = data.frame(tax_metatrans_summary, method=rep("metatranscriptomics",8),metadata[,c("subject_id","fiber")] )
tax_16S_summary       = t(prop.table(apply(tax_16S[,2:9] , 2, tapply, tax_16S$manual_tax, sum),2))
tax_16S_summary       = data.frame(tax_16S_summary , Archaea = rep(0,8), method=rep("16S rRNA genes seq.",8), metadata[,c("subject_id","fiber")] )

rownames(tax_16S_summary) = rownames(tax_metatrans_summary) = colnames(SubCatKegg)

metatrans_alimintest = list(
	tax_metatrans_summary = tax_metatrans_summary, 
	tax_16S_summary       = tax_16S_summary, 
	SubCatKegg            = SubCatKegg, 
	metadata              = metadata
	)

save(metatrans_alimintest, file = "data/metatrans_alimintest.RData")
