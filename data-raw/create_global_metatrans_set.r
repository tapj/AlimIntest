

# prepare data
ind= c("79","96","99","100","126","151","175")

## read blast, kegg linked with sample id and reads
blast_all = NULL
kegg_all = NULL

for (i in ind) {


blast = read.table(paste0(i,"f_V3_clean.blastx"), sep="\t")
kegg = read.table(paste0(i,"f_V3_clean-U-Ecological_Table"), sep="\t")

blast = cbind(rep(i, dim(blast)[1]) , blast)
kegg = cbind(rep(i, dim(kegg)[1]) , kegg)

blast_all = rbind(blast_all, blast)
kegg_all = rbind(kegg_all, kegg)


}



### remove spurious hit adn clean table

blast_all = blast_all[blast_all$V11 < 0.00001,]

blast_all_clean = cbind(blast_all[1:3], 
  do.call(
      rbind,
      lapply(
        strsplit(
          as.character(blast_all$V1), split="_"),
          function(x){ 
                x[1] 
              }
            )
          )
         )
         
         
         
names(blast_all_clean) =c("ind","read.trim","ncbi_id","raw_read_id")



kegg_all_clean = cbind(kegg_all, 
  do.call(
      rbind,
      lapply(
        strsplit(
          as.character(kegg_all$V1), split="_"),
          function(x){ 
                x[1] 
              }
            )
          )
         )
names(kegg_all_clean)[c(1,9)] = c("ind_id","raw_read_id")


#### read KO sub categories correspondence
Kegg_subcat_corresp = read.table("Kegg_subcat_corresp.txt", sep=" ", colClasses = c("character","character"))
names(Kegg_subcat_corresp) = c("SubCat","KO")

#kegg_all_clean = merge(kegg_all_clean, Kegg_subcat_corresp, by.x="V3", by.y="KO", all=TRUE)

ncbi_kegg_clean = merge(blast_all_clean, kegg_all_clean, all=TRUE)


ncbi_kegg_clean = cbind(ncbi_kegg_clean, 
  do.call(
      rbind,
      lapply(
        strsplit(
          as.character(ncbi_kegg_clean$ncbi_id), split="\\|"),
          function(x){ 
                x[4] 
              }
            )
          )
         )

names(ncbi_kegg_clean)[13] = "accnum"



## read NCBI prot id linked to taxonomy



ncbi_taxo      = read.table("prot_taxo.txt", sep="\t")
ncbi_kegg_taxo = merge(ncbi_kegg_clean, ncbi_taxo, by.x="accnum", by.y="V1", all.x=TRUE)



## read cazy linked with contig and singlet
cazy = read.csv2("metatrans_CAZY.csv")


## read reads linked to contig
read_to_contig = read.table("pool_tag_outfile1.txt", sep="\t")
read_to_contig = cbind(read_to_contig[c(1,3)], 
  do.call(
      rbind,
      lapply(
        strsplit(
          as.character(read_to_contig$V1), split="_"),
          function(x){ 
                rbind(c(x[1],x[2])) 
              }
            )
          )
         )
         
         
names(read_to_contig) = c("read_id", "contig_id", "ind","raw_read_id")

read_to_contig$ind = gsub(" ", "",(read_to_contig$ind))

#### merge with singlets

singlets = readLines("Alimintest.singlets.cap.fasta")
singlets = singlets[grep(">",singlets)]
singlets = gsub(">","", singlets)
singlets = gsub(" ", "", singlets)

singlets = cbind(singlets, singlets,
  do.call(
      rbind,
      lapply(
        strsplit(
          as.character(singlets), split="_"),
          function(x){ 
                rbind(c(x[1],x[2])) 
              }
            )
          )
         )
         
colnames(singlets) = c("read_id", "contig_id", "ind","raw_read_id")

read_to_contig = rbind(read_to_contig, as.data.frame(singlets))



# link all data into one catalog

read.cazy = merge(read_to_contig, cazy, by.x="contig_id", by.y="CONTIG", all.x = TRUE)

read.cazy.kegg.ncbi = merge(read.cazy, ncbi_kegg_taxo, by="raw_read_id", all=TRUE); gc()
read.cazy.kegg.ncbi = unique(read.cazy.kegg.ncbi)

Eukaryota_idx = grep("Eukaryota", read.cazy.kegg.ncbi$V2.y)

read.cazy.kegg.ncbi = read.cazy.kegg.ncbi[-Eukaryota_idx,] #remove Eukaryota

dim(read.cazy.kegg.ncbi)

write.table(read.cazy.kegg.ncbi, file="metatrans_read.cazy.kegg.ncbi.txt", sep="\t")

#Analysis
library(ade4)


# ncbi_cazy_kegg_table = table(as.character(read.cazy.kegg.ncbi$V2.y), paste(read.cazy.kegg.ncbi$Functional.annotation, read.cazy.kegg.ncbi$V6))

# ncbi_cazy_table = as.data.frame.matrix(table(read.cazy.kegg.ncbi$V2.y, read.cazy.kegg.ncbi$Functional.annotation))

# ncbi_kegg_table = as.data.frame.matrix(table(read.cazy.kegg.ncbi$V2.y, read.cazy.kegg.ncbi$V6))

# ncbi_cazy_kegg_table = cbind(ncbi_cazy_table,ncbi_kegg_table)
# ncbi_cazy_kegg_table = ncbi_cazy_kegg_table[apply(ncbi_cazy_kegg_table , 1, sum) > 0,apply(ncbi_cazy_kegg_table , 2, sum) > 0]
# dim(ncbi_cazy_kegg_table)

# ncbi_cazy_kegg_table = as.data.frame.matrix(ncbi_cazy_kegg_table)

# meta.coa = dudi.coa(as.data.frame(t(ncbi_cazy_kegg_table)) , scannf=F, nf=3)




### Co-inerta subcat KEGG and ncbi taxonomy


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

metatranscriptomics_alimintest = list(tax_metatrans_summary, tax_16S_summary, SubCatKegg, metadata)

save(metatranscriptomics_alimintest, file = "data/metatranscriptomics_alimintest.RData")

# analysis

library(ggplot2)
library(ade4)


tax_levels = c(

    "Bacteroides",
    "O_Bacteroidetes",
    "Bifidobacterium",
    "O_Actinobacteria",
    "Clostridiaceae",
    "Eubacteriaceae",
    "Lachnospiraceae",
    "Ruminococcaceae",
    "O_Firmicutes",
    "O_Proteobacteria",
    "O_Bacteria",
    "Archaea"

) #useful to reorder the taxonomic level for the bar chart plot


metatrans = rbind(

	tax_metatrans_summary[,c(tax_levels, "method","subject_id", "fiber")],
	tax_16S_summary[,c(tax_levels, "method","subject_id", "fiber") ]

)
 
 
metatrans = melt(metatrans, id.vars=c("method","subject_id","fiber"))
 
taxcolors = c("blue","deepskyblue1","yellow","yellow3", "hotpink","hotpink1","hotpink2","hotpink3",
"hotpink4", "palegreen", "grey","black")
 
#pdf("~/storage/Dropbox/AlimIntest_reloaded/Article/revised_manuscript/Supplementary/FigS6_metatrans_tax.pdf", h=7, w=7)
 
ggplot(metatrans) + geom_bar(aes(y=value, x=fiber, fill=variable), stat="identity") + facet_grid(method~subject_id) + scale_x_discrete(label=c("10g","40g")) + scale_fill_manual("Taxonomy",values=taxcolors) + xlab("fiber amount per day") + ylab("proportions")
 
#dev.off()
 
tax        = t(tax_metatrans_summary[1:12])
#tax        = t(tax_16S_summary[1:11])




# tax_subCatKegg.coi = coinertia(
# dudi.pca(as.data.frame(t(log10(tax+0.00001))),scannf=F, nf=3),
# dudi.pca(as.data.frame(t(log10(SubCatKegg[1:15,]+0.00001))),scannf=F, nf=3),
# scannf=F, nf=2)

# plot(tax_subCatKegg.coi)

# tax_subCatKegg.coi.bca = coinertia(
# bca(dudi.pca(as.data.frame(t(log10(tax+0.00001))),scannf=F, nf=3),fac=metadata$fiber, scannf=F, nf=1),
# bca(dudi.pca(as.data.frame(t(log10(SubCatKegg[1:15,]+0.00001))),scannf=F, nf=3),fac=metadata$fiber, scannf=F, nf=1),
# scannf=F, nf=2)



# plot(tax_subCatKegg.coi$co)
# text(tax_subCatKegg.coi$li,label=row.names(tax_subCatKegg.coi$li))
# text(tax_subCatKegg.coi$co,label=row.names(tax_subCatKegg.coi$co), col="red")



wit1 <- wca(dudi.pca(as.data.frame(t(log10(tax+0.00001))),scannf=F, nf=3), metadata$fiber, scan = FALSE, scal = "total")
wit2 <- wca(dudi.pca(as.data.frame(t(log10(SubCatKegg[1:15,]+0.00001))),scannf=F, nf=3), metadata$fiber, scan = FALSE, nf = 2)
kta1 <- ktab.within(wit1, colnames = metadata$sample_id2)
kta2 <- ktab.within(wit2, colnames = metadata$sample_id2)
statico1 <- statico(kta1, kta2, scan = FALSE)

plot(statico1)






