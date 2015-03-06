# Generate alimintestData from raw data

metadata = read.csv2("data-raw/metadata.csv", header=TRUE, row.names=1)
otu      = read.csv2("data-raw/otu.csv", header=TRUE, row.names=1, check.names=FALSE)
tax      = read.csv2("data-raw/tax.csv", header=TRUE, row.names=1)
basalqpcr = read.csv2("data-raw/metadata_basal_qpcr.csv", header=TRUE, row.names=1)

alimintestData=list(metadata=metadata, otu=otu, tax=tax, basalqpcr=basalqpcr)

save(alimintestData, file="data/alimintestData.RData")


# prepare data for replicate analysis



sanger_454   = read.csv2("data-raw/sanger.454.csv",   header = TRUE, row.names = 1)
replicat_454 = read.csv2("data-raw/replicat.454.csv", header = TRUE, row.names = 1)

sanger_454$tax       = as.character(sanger_454$tax)
replicat_454$tax     = as.character(replicat_454$tax)
replicat_454$conf    = as.numeric(as.character(replicat_454$conf))
replicat_454$samples = as.character(replicat_454$samples)

save(sanger_454, file="data/sanger_454.RData")
save(replicat_454, file="data/replicat_454.RData")


######################
