### This code is meant to create a large 'master' table of "project" somatic WGS data
### using SNV and INDEL vcf data, all mutations from all samples will be combined

## Import libraries ------------------------------------------------------------
library(openxlsx)
library(dplyr)

## Data wrangling --------------------------------------------------------------
# Import merged vcf data
# x <- importVCFfiles(vcfFiles = '/path/to/merged_annotated.vcf')

# Read annotated excel file
data <- read.delim('/path/to/mutations_annotation.tsv', header = TRUE, sep = '\t')

# Replace "None" to make data easier to read
data[data == "None"] <- "."

# Re-order columns and ignore duplicate germline columns
# new_order <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", SAMPLE ID COLUMNS --REDACTED)

# Re-order annotated columns (previous version)
# new_order <- c("Sample_ID", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "SYMBOL","Consequence",
  #               "Sample_GT", "Sample_DP", "Sample_AAF", "Sample_AD1", "Sample_AD2",
  #              "AC", "AF","IMPACT","gnomAD_exomes_AF", "gnomAD_genomes_AF", "gnomAD_genomes_NFE_AF",
  #              "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG",
  #             "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
  #              "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_SYMBOL",
  #              "MES_NCSS_downstream_acceptor", 
  #              "MES_NCSS_downstream_donor", "MES_NCSS_upstream_acceptor",
  #              "MES_NCSS_upstream_donor", "MES_SWA_acceptor_alt", "MES_SWA_acceptor_diff", 
  #              "MES_SWA_acceptor_ref", "MES_SWA_acceptor_ref_comp", "MES_SWA_donor_alt",
  #              "MES_SWA_donor_diff", "MES_SWA_donor_ref", "MES_SWA_donor_ref_comp",
  #              "ALFA_European_AF", "AN", "Allele", "Amino_acids",
  #              "BIOTYPE", "BayesDel_addAF_pred", "BayesDel_addAF_rankscore", 
  #              "BayesDel_addAF_score", "BayesDel_noAF_pred", "BayesDel_noAF_rankscore",
  #              "BayesDel_noAF_score", "CAVA_ALTANN", "CAVA_ALTCLASS", "CAVA_ALTFLAG", 
  #              "CAVA_ALTSO", "CAVA_CLASS", "CAVA_CSN", "CAVA_C", "CAVA_P", "CAVA_GENE",
  #              "CAVA_GENEID", "CAVA_IMPACT", "CAVA_LOC", "CAVA_PROTALT", "CAVA_PROTPOS",
  #              "CAVA_PROTREF", "CAVA_SO", "CAVA_TRANSCRIPT", "CAVA_TRINFO", "CAVA_TYPE",
  #              "CDS_position", "ClinVar", "ClinVar_CLNDISDB", "ClinVar_CLNDN", 
  #              "ClinVar_CLNREVSTAT", "ClinVar_CLNSIG", "ClinVar_CLNSIGCONF", 
  #              "ClinVar_CLNSIGINCL", "ClinVar_CLNVC", "ClinVar_CLNVI", "Codons",
  #              "DISTANCE", "DP", "ECNT", "ENSP", "EXON", "ExAC_Adj_AF",
  #              "Existing_variation", "FLAGS", "Feature", "Feature_type", "Gene", "HGNC_ID",
  #              "HGVS_OFFSET", "HGVSc", "HGVSp", "INTRON", "IN_PON", "LOEUF", 
  #              "MANE_PLUS_CLINICAL", "MANE_SELECT",
  #              "MaveDB_nt", "MaveDB_pro", "MaveDB_score", "MaveDB_urn", "MaxEntScan_alt",
  #              "MaxEntScan_diff", "MaxEntScan_ref", "NLOD", "NMD", "N_ART_LOD", "POP_AF",
  #              "P_CONTAM", "P_GERMLINE", "PrimateAI", "Protein_position", "REVEL_rankscore", 
  #              "REVEL_score", "RPA", "RU", "SOURCE", "STR", "STRAND", "SYMBOL_SOURCE", 
  #              "TLOD",
  #              "UK10K_AF", "VARIANT_CLASS", "am_class", "am_pathogenicity", "cDNA_position",
  #              "pLI_gene_value") 

data$Sample_AAF <- as.numeric(data$Sample_AAF)
summary(data$Sample_AAF) # Check for NA values

data[is.na(data$Sample_AAF), ]

data <- data[!is.na(data$Sample_AAF), ]

data <- data %>% mutate(Sample_VAF = 1 - Sample_AAF)

# Re-order annotated columns
new_order <- c("Sample_ID", "CHROM", "POS", "REF", "ALT", "SYMBOL","Consequence",
              "Sample_GT", "Sample_DP", "Sample_VAF", "Sample_AD1", "Sample_AD2",
             "AC","IMPACT","AN", "Allele", "Amino_acids","CAVA_ALTANN", 
             "CAVA_CSN", "CAVA_C", "CAVA_P", "CAVA_GENE",
             "CAVA_GENEID", "CAVA_IMPACT", "CAVA_LOC", "CAVA_SO", "CAVA_TRANSCRIPT","CAVA_TYPE",
             "gnomAD_exomes_AF", "gnomAD_genomes_AF", "gnomAD_genomes_NFE_AF",
             "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG",
             "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
             "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_SYMBOL",
             "ALFA_European_AF", 
             "BIOTYPE", "BayesDel_addAF_pred", "BayesDel_addAF_rankscore", 
             "BayesDel_addAF_score", "BayesDel_noAF_pred", "BayesDel_noAF_rankscore",
             "BayesDel_noAF_score", "CAVA_ALTCLASS", "CAVA_ALTFLAG", 
             "CAVA_ALTSO", "CAVA_CLASS", "CAVA_PROTALT", "CAVA_PROTPOS",
             "CAVA_PROTREF", "CAVA_TRINFO",
             "CDS_position", "ClinVar", "ClinVar_CLNDISDB", "ClinVar_CLNDN", 
             "ClinVar_CLNREVSTAT", "ClinVar_CLNSIG", "ClinVar_CLNSIGCONF", 
             "ClinVar_CLNSIGINCL", "ClinVar_CLNVC", "ClinVar_CLNVI", "Codons",
             "DISTANCE", "DP", "ECNT", "ENSP", "EXON", "ExAC_Adj_AF",
             "Existing_variation", "FLAGS", "Feature", "Feature_type", "Gene", "HGNC_ID",
             "HGVS_OFFSET", "HGVSc", "HGVSp", "INTRON", "IN_PON", "LOEUF", 
             "MANE_PLUS_CLINICAL", "MANE_SELECT",
             "MaveDB_nt", "MaveDB_pro", "MaveDB_score", "MaveDB_urn", "MaxEntScan_alt",
             "MaxEntScan_diff", "MaxEntScan_ref", "NLOD", "NMD", "N_ART_LOD", "POP_AF",
             "P_CONTAM", "P_GERMLINE", "PrimateAI", "Protein_position", "REVEL_rankscore", 
             "REVEL_score", "RPA", "RU", "SOURCE", "STR", "STRAND", "SYMBOL_SOURCE", 
             "TLOD",
             "UK10K_AF", "VARIANT_CLASS", "am_class", "am_pathogenicity", "cDNA_position",
             "pLI_gene_value")

rearrangedData <- data %>% select(all_of(new_order))

## Write rearrangedData to an Excel file ---------------------------------------
write.xlsx(rearrangedData, file = "/path/to/output.xlsx")
