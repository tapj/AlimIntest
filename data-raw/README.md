# _AlimIntest_ raw dataset

## Metadata
	
    * idx:	INDEX
    * Sample_ID; Sample_ID2:	Sample identifications
    * Time point:	Time point corresponding to study design illustrated on Fig. S1
    * Subject_Id:	Subject idenfication
    * Diet_run:	Diet run "r10_40" correspond to 10g fiber per day between time point n°2 and time point n°3 and 40g fiber per day between time point n°4 and time point n°5; "r40_10" correspond to the cross over to "r10_40". See Fig. S1
    * Seq_filtered:	Sequence after quality checking from LOTUS pipeline
    * acetate; propionate; isobutyrate; butyrate; isovalerate; valerate; isocaproate; caproate:	Short Chain Fatty Acid measurement  (mM) by gas-liquid chromatography
    * log_all_Bacteria; log_Cleptum_group; log_Bcoccoides_group; log_Bacteroides_Prevotella; log_Bifidobacteria; log_Ecoli; log_Lactobacillus:	Quantitative PCR assays expressed in Log10 of number of equivalent bacterial cells
    * fecal_water_pH:	Fecal water pH measurement
    * Comet_assay:	Single Cell Gel Electrophoresis assay expressed in % of DNA in the tail
		
## otu
	
	* OTUs table from LOTUS pipeline (clustered at 97% similarity), column names correspond to Sample_ID, row names correpond to OTU id	
		
## tax	

	* Taxonomic identification for each OTU, taxonomy uncertainties are illustrated by "?"	
		
## sequences\_files_location

	* 454 Sequences files location for 16S rRNA gene amplicon and metatranscriptomic assays with corresponding Sample_ID and accession number	
