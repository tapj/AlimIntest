#' AlimIntest data
#'
#' This dataset contains all metadata and microbial data from the 
#' AlimIntest project.
#'
#'
#' @section Variables:
#'
#' \itemize{
#' \item \code{metadata}: study design with SCFA, qPCR and comet assays
#'			\itemize{
#'    \item idx:	INDEX
#'    \item Sample_ID; Sample_ID2:	Sample identifications
#'    \item Time point:	Time point corresponding to study design illustrated on Fig. S1
#'    \item Subject_Id:	Subject identification
#'    \item Diet_run:	Diet run "r10_40" correspond to 10g fiber per day between time point n°2 and time point n°3 and 40g fiber per day between time point n°4 and time point n°5; "r40_10" correspond to the cross over to "r10_40". See Fig. S1
#'    \item Seq_filtered:	Sequence after quality checking from LOTUS pipeline
#'    \item acetate; propionate; isobutyrate; butyrate; isovalerate; valerate; isocaproate; caproate:	Short Chain Fatty Acid measurement  (mM) by gas-liquid chromatography
#'    \item log_all_Bacteria; log_Cleptum_group; log_Bcoccoides_group; log_Bacteroides_Prevotella; log_Bifidobacteria; log_Ecoli; log_Lactobacillus:	Quantitative PCR assays expressed in Log10 of number of equivalent bacterial cells
#'    \item fecal_water_pH:	Fecal water pH measurement
#'    \item Comet_assay:	Single Cell Gel Electrophoresis assay expressed in % of DNA in the tail
#'			}
#' \item \code{otu}: Microbial OTU table
#' \item \code{tax}: OTU Greengenes taxonomy
#' \item \code{basalqpcr}: extra qpcr data which correspond to time point 1 and time point 6 (see Fig. S1)
#' }
#' @docType data
#' @name alimintestData
#' @usage alimintestData
#' @format A list of 4 data frames
#' @examples
#' head(alimintestData$metadata)
NULL

#' AlimIntest Sanger VS 454 sequencing data
#'
#' This dataset contains  Sanger VS 454 sequencing data from the same DNA sample
#' AlimIntest project.
#'
#'
#' @section Variables:
#'
#' \itemize{
#' \item \code{samples}: sample replicates id
#' \item \code{technology}: sequencing technology
#' \item \code{run}: how replicates were sequenced
#' \item \code{tax}: taxonomy assigned by RDP classifier at genus levels
#' \item \code{conf}: confidence score from the RDP classifier. 
#' }
#' @docType data
#' @name sanger_454
#' @usage sanger_454
#' @format 1 data frame with 15504 rows (reads) and 5 columns (variables)
#' @examples
#' head(sanger_454)
NULL


#' 454 sequencing data with 4 replicate samples from two manufacturers
#'
#' This dataset contains  Sanger VS 454 sequencing data from the same DNA sample
#' AlimIntest project.
#'
#'
#' @section Variables:
#'
#' \itemize{
#' \item \code{samples}: sample id
#' \item \code{manufacturer}: sequencing 454 manufacturer
#' \item \code{tax}: taxonomy assigned by RDP classifier at genus levels
#' \item \code{conf}: confidence score from the RDP classifier. 
#' }
#' @docType data
#' @name replicat_454
#' @usage replicat_454
#' @format 1 data frame with 87782 rows (reads) and 4 columns (variables)
#' @examples
#' head(replicat_454)
NULL

#' Metatranscriptomics AlimIntest data
#'
#' This dataset contains Metatranscriptomics data from the 
#' AlimIntest project.
#'
#'
#' @section Variables:
#'
#' \itemize{
#' \item \code{tax_metatrans_summary}: taxonomic table from metatranscriptomics sequencing
#' \item \code{tax_16S_summary}: taxonomic table from 16S rRNA genes reads
#' \item \code{SubCatKegg}: KEGG categories from metatranscriptomics sequencing
#' \item \code{metadata}: minimal set of metadata for each sample analysed by metatranscriptomics
#' }
#' @docType data
#' @name metatrans_alimintest
#' @usage metatrans_alimintest
#' @format A list of 4 data frames
#' @examples
#' head(metatrans_alimintest$metadata)
NULL

#' Basal diet AlimIntest data
#'
#' This dataset contains basal diet from individual before and after 
#' inclusion into the AlimIntest study
#'
#'
#' @section Variables:
#'
#' \itemize{
#' \item \code{subject}: Subject identification
#' \item \code{visit}: before or after nutritional intervention
#' \item \code{day}: nb of days before or after nutritional intervention
#' \item \code{food}: 124 food items, "1" means food item included into diet by the subject
#' }
#' @docType data
#' @name basal_diet
#' @usage basal_diet
#' @format a data frame with 114 rows and 127 columns
#' @examples
#' head(basal_diet)
NULL


#' calculate delta for a variable
#'
#' This functions permit to calculate difference before and after diet 
#'
#'
#' @param Subject_Id Subject_Id should be a factor
#' @param time.point time.point should be a factor with levels 2,3,4,5
#' @param Diet_run should be a factor with levels r40-10 and r10-40
#' @param metadata data to be processed
#' @param var_name the variable name to be processed inside metadata
#' @return res result
#' @examples
#' data(alimintestData)
#' calculate_delta(Subject_Id=alimintestData$metadata$Subject_Id,time.point=alimintestData$metadata$Time.point,
#'                Diet_run=alimintestData$metadata$Diet_run,metadata=alimintestData$metadata,var_name="acetate")
#' @export

calculate_delta=function(Subject_Id,time.point,Diet_run,metadata,var_name="variable") {

res=NULL

variable=metadata[,var_name]

for(s in levels(Subject_Id)) {

	if(Diet_run[Subject_Id==s][1]=="r40-10") {
		
		before.diet.40 = variable[which(Subject_Id==s & time.point=="2")]
		before.diet.10 = variable[which(Subject_Id==s & time.point=="4")]
		after.diet.40  = variable[which(Subject_Id==s & time.point=="3")]
		after.diet.10  = variable[which(Subject_Id==s & time.point=="5")]
		
	} else {
		
		before.diet.10 = variable[which(Subject_Id==s & time.point=="2")]
		before.diet.40 = variable[which(Subject_Id==s & time.point=="4")]
		after.diet.10  = variable[which(Subject_Id==s & time.point=="3")]
		after.diet.40  = variable[which(Subject_Id==s & time.point=="5")]
	
	}

delta.10=before.diet.10-after.diet.10
delta.40=before.diet.40-after.diet.40

dd=data.frame(Subject_Id=rep(s,2), fiber.diet=c("delta.10","delta.40"), variable=rep(var_name,2), delta=c(delta.10, delta.40))

res=rbind(res,dd)



}

return(res)
}

