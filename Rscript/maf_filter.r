# Rscripts maf_filter.r inputfile outputfile1 outputfile2 outputfile3 outputfile4 genefile
library(dplyr)
library(stringr)

inputfile = commandArgs(trailingOnly=TRUE)[1]
outputfile1 = commandArgs(trailingOnly=TRUE)[2]
outputfile2 = commandArgs(trailingOnly=TRUE)[3]
outputfile3 = commandArgs(trailingOnly=TRUE)[4]
input_intervar = commandArgs(trailingOnly=TRUE)[5]
input_intervar_raw = commandArgs(trailingOnly=TRUE)[6]

#genes <- read.csv(genefile, sep="\t")
#onco_genes <- genes$Hugo.Symbol

# filtering variant and add vaf
# input : _vep_oncoanno_raw.maf
matrix <- read.csv(inputfile, header=T, sep="\t")
intervar <- read.csv(input_intervar, header = T, sep="\t")
intervar_raw <- read.csv(input_intervar_raw, header = T, sep="\t")
#merged <- merge(matrix, intervar[c(2:6,13:14)], by.x=c('Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Hugo_Symbol'), by.y=c('Start', 'End', 'Ref', 'Alt', 'Ref.Gene'))
#merged <- merge(matrix, intervar[c(2:9,13:14)], by.x=c('Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Hugo_Symbol'), by.y=c('Start', 'End', 'Ref', 'Alt', 'Gene.ensGene'))
#merged <- merge(matrix, intervar[c(2:9,13:14)], by.x=c('Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'), by.y=c('Start', 'End', 'Ref', 'Alt'), all.x = TRUE)

matrix <- matrix %>%
  mutate(
    Gene_vcfid = sub(":.*$", "", vcf_id)
  )

matrix <- matrix %>%
  mutate(
    CLNREVSTAT = case_when(
      is.na(CLNREVSTAT) | CLNREVSTAT == "" ~ CLNREVSTAT,
      CLNREVSTAT == "practice_guideline" ~ "four",
      CLNREVSTAT == "reviewed_by_expert_panel" ~ "three",
      CLNREVSTAT == "criteria_provided,_multiple_submitters,_no_conflicts" ~ "two",
      CLNREVSTAT == "criteria_provided,_conflicting_classifications" ~ "one",
      CLNREVSTAT == "criteria_provided,_single_submitter" ~ "one",
      TRUE ~ "zero"
    )
  )

intervar <- intervar %>%
  mutate(
    Gene_otherinfo = sub(":.*$", "", Otherinfo)
  )
 
intervar_raw <- intervar_raw %>%
  mutate(
    Gene_otherinfo3 = sub(":.*$", "", Otherinfo3)
  )
 
 cols_keep <- c(
  "Chr", "Start", "End", "Ref", "Alt",
  "SIFT_score",
  "Polyphen2_HDIV_score",
  "Polyphen2_HVAR_score",
  "MutationTaster_score",
  "MutationAssessor_score",
  "FATHMM_score",
  "PROVEAN_score",
  "VEST4_score",
  "MetaLR_score",
  "MetaRNN_score",
  "M-CAP_score",
  "REVEL_score",
  "MutPred_score",
  "MVP_score",
  "MPC_score",
  "PrimateAI_score",
  "DEOGEN2_score",
  "BayesDel_addAF_score",
  "BayesDel_noAF_score",
  "ClinPred_score",
  "DANN_score",
  "fathmm-MKL_coding_score",
  "Gene_otherinfo3"
)

intervar_raw_sub <- intervar_raw %>%
  select(any_of(cols_keep))

intervar_key1 <- intervar %>%
  select(X.Chr, Start, End, Ref, Alt, Gene_otherinfo, everything()) %>%
  distinct(X.Chr, Start, End, Ref, Alt, Gene_otherinfo, .keep_all = TRUE)

intervar_merged <- intervar_key1 %>%
  left_join(
    intervar_raw_sub,
    by = c(
      "Start" = "Start",
      "End"   = "End",
      "Ref"   = "Ref",
      "Alt"   = "Alt",
	  "Gene_otherinfo" = "Gene_otherinfo3"
    )
  )

merged <- matrix %>%
  left_join(
    intervar_merged,
    by = c(
      'Start_Position'    = 'Start',
      #'End_Position'      = 'End',
      'Reference_Allele'  = 'Ref',
      'Tumor_Seq_Allele2' = 'Alt',
	  "Gene_vcfid" = "Gene_otherinfo"
    )
  )



# numerical conversion
merged$ToMMo_AF <- as.numeric(merged$ToMMo_AF)

# variant classification
matrix1 <- merged %>%
  mutate(
    # ========= flag =========
    rev_good = CLNREVSTAT %in% c("one", "two", "three", "four"),
    af_low   = is.na(ToMMo_AF) | ToMMo_AF < 0.001,
    af_high  = !is.na(ToMMo_AF) & ToMMo_AF >= 0.001,

    onc_pos  = str_detect(ONCOGENIC, "Oncogenic|Likely Oncogenic"),
    onc_lpos = str_detect(ONCOGENIC, "Likely Oncogenic"),
    onc_neu  = str_detect(ONCOGENIC, "Neutral|Likely Neutral|Inconclusive"),

    ivar_path   = str_detect(InterVar..InterVar.and.Evidence,
                             "Pathogenic|Likely pathogenic"),
    ivar_benign = str_detect(InterVar..InterVar.and.Evidence,
                             "Benign|Likely benign")
  ) %>%
  mutate(
    pipeline_class = case_when(

      # ======================================================
      # Pathogenic CLNSIG
      # ======================================================
      CLNSIG %in% c("Pathogenic","Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic","Pathogenic|other") &
        rev_good ~ "Pathogenic/Likely_pathogenic",

      CLNSIG %in% c("Pathogenic","Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic","Pathogenic|other") &
        CLNREVSTAT == "zero" & af_low &
        (onc_pos | ivar_path) ~ "Pathogenic/Likely_pathogenic",

      CLNSIG %in% c("Pathogenic","Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic","Pathogenic|other") &
        CLNREVSTAT == "zero" & af_low ~ "VUS",

      CLNSIG %in% c("Pathogenic","Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic","Pathogenic|other") &
        af_high ~ "VUS",

      # ======================================================
      # drug_response / risk_factor / other / not_provided
      # ======================================================
      # CLNSIG %in% c("drug_response","risk_factor","other","not_provided") &
      #   (onc_pos | onc_lpos) ~ "Pathogenic/Likely_pathogenic",
      CLNSIG %in% c("drug_response","risk_factor","other","not_provided", "Uncertain_significance|drug_response", "Uncertain_significance/Uncertain_risk_allele") &
        (onc_pos | ivar_path) & af_low ~ "Pathogenic/Likely_pathogenic",

      CLNSIG %in% c("drug_response","risk_factor","other","not_provided", "Uncertain_significance|drug_response", "Uncertain_significance/Uncertain_risk_allele") &
        (onc_pos | ivar_path) & af_high ~ "VUS",

      CLNSIG %in% c("drug_response","risk_factor","other","not_provided", "Uncertain_significance|drug_response", "Uncertain_significance/Uncertain_risk_allele") &
        (onc_neu | ivar_benign) ~ "Benign/Likely_benign",

      # CLNSIG %in% c("drug_response","risk_factor","other","not_provided") &
      #   ONCOGENIC == "Unknown" ~ "VUS",

      CLNSIG %in% c("drug_response","risk_factor","other","not_provided", "Uncertain_significance|drug_response", "Uncertain_significance/Uncertain_risk_allele") ~ "VUS",

      # ======================================================
      # Benign CLNSIG
      # ======================================================
      CLNSIG %in% c("Benign","Likely_benign","Benign/Likely_benign") &
        rev_good ~ "Benign/Likely_benign",

      CLNSIG %in% c("Benign","Likely_benign","Benign/Likely_benign") &
        CLNREVSTAT == "zero" & af_high ~ "Benign/Likely_benign",

      CLNSIG %in% c("Benign","Likely_benign","Benign/Likely_benign") &
        CLNREVSTAT == "zero" & af_low &
        (onc_neu | ivar_benign) ~ "Benign/Likely_benign",

      CLNSIG %in% c("Benign","Likely_benign","Benign/Likely_benign") &
        CLNREVSTAT == "zero" ~ "VUS",

      # ======================================================
      # Uncertain significance complex
      # ======================================================
      # CLNSIG %in% c("Uncertain_significance|drug_response",
      #               "Uncertain_significance/Uncertain_risk_allele") &
      #   (onc_pos | ivar_path) ~ "Pathogenic/Likely_pathogenic",
      #
      # CLNSIG %in% c("Uncertain_significance|drug_response",
      #               "Uncertain_significance/Uncertain_risk_allele") ~ "VUS",

      # ======================================================
      # Uncertain_significance solo
      # ======================================================
      CLNSIG == "Uncertain_significance" & rev_good ~ "VUS",

      CLNSIG == "Uncertain_significance" &
        CLNREVSTAT == "zero" & af_low &
        (onc_pos | ivar_path) ~ "Pathogenic/Likely_pathogenic",

      CLNSIG == "Uncertain_significance" ~ "VUS",

      # ======================================================
      # Conflicting classifications
      # ======================================================
      CLNSIG == "Conflicting_classifications_of_pathogenicity" &
        str_detect(CLNSIGCONF, "Pathogenic|Likely_pathogenic") &
        str_detect(CLNSIGCONF, "Benign|Likely_benign") ~ "VUS",

      CLNSIG == "Conflicting_classifications_of_pathogenicity" &
        str_detect(CLNSIGCONF, "Pathogenic|Likely_pathogenic") &
        (onc_pos | ivar_path) ~ "Pathogenic/Likely_pathogenic",

      CLNSIG == "Conflicting_classifications_of_pathogenicity" &
        str_detect(CLNSIGCONF, "Pathogenic|Likely_pathogenic") &
        !(onc_pos | ivar_path) ~ "VUS",

      CLNSIG == "Conflicting_classifications_of_pathogenicity" &
        str_detect(CLNSIGCONF, "Benign|Likely_benign") &
        (onc_pos | ivar_path) ~ "VUS",

      CLNSIG == "Conflicting_classifications_of_pathogenicity" &
        str_detect(CLNSIGCONF, "Benign|Likely_benign") ~ "Benign/Likely_benign",

      # ======================================================
      # CLNSIG blank
      # ======================================================
      (is.na(CLNSIG) | CLNSIG == "") &
        !(onc_pos | ivar_path) ~ "VUS",

      (is.na(CLNSIG) | CLNSIG == "") &
        (onc_pos | ivar_path) & af_low ~ "Pathogenic/Likely_pathogenic",

      (is.na(CLNSIG) | CLNSIG == "") &
        (onc_pos | ivar_path) & af_high ~ "VUS",

      # ======================================================
      # fallback
      # ======================================================
      TRUE ~ "undefined"
    )
  )

# VUS prioritization
matrix2 <- matrix1 %>%
  mutate(
    # ===== numeric conversion =====
    CADD_phred              = as.numeric(CADD_phred),
    Polyphen2_HDIV_score    = as.numeric(Polyphen2_HDIV_score),
    REVEL_score             = as.numeric(REVEL_score),
    BayesDel_noAF_score     = as.numeric(BayesDel_noAF_score),
    MetaSVM_score           = as.numeric(MetaSVM_score),
    PrimateAI_score         = as.numeric(PrimateAI_score)
  ) %>%
  mutate(
    # ===== pipeline_class_vus_prior =====
    pipeline_class_vus_prior = case_when(
      pipeline_class %in% c("Pathogenic/Likely_pathogenic",
                            "Benign/Likely_benign") ~ pipeline_class,

      pipeline_class == "VUS" & (
        CADD_phred >= 28.1 |
        #Polyphen2_HDIV_score >= 0.999 |
        REVEL_score >= 0.773 |
        BayesDel_noAF_score >= 0.27 |
        MetaSVM_score >= 0.8 |
        PrimateAI_score >= 0.867 |
        SpliceAI_cutoff == "PASS"
      ) ~ "VUS_high",

      pipeline_class == "VUS" ~ "VUS",

      TRUE ~ NA_character_
    )
  )


# outfile1
write.table(matrix1, outputfile1, na="", row.names=F, quote=F, sep="\t")

# outfile2
write.table(matrix2, outputfile2, na="", row.names=F, quote=F, sep="\t")

# outfile3
write.table(merged, outputfile3, na="", row.names=F, quote=F, sep="\t")
