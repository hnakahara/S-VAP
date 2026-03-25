container: "docker://condaforge/mambaforge:24.9.2-0"

# true or false
GRCh38 = "false"

filetype = "vcf"

if GRCh38 == "false":
    print("Genome version: GRCh37")
elif GRCh38 == "true":
    print("Genome version: GRCh38")
else:
    print("GRCh38 is unknown parameter")

################################
# Rule all
################################

rule all:
    input:
        "results/libcrypto.so.1.0.0_copy.done",
        "results/11final-output/merge_variant_class_vus_prior.maf"

################################
# Modules
################################

include: "rules/preprocessing.smk"
include: "rules/annotation_from_vcf.smk"

# MAF filter
if GRCh38 == "false":
    include: "rules/maf_filter_hg19_from_vcf.smk"
else:
    include: "rules/maf_filter_hg38_from_vcf.smk"