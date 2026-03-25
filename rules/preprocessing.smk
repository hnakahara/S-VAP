import pandas as pd
import numpy as np
import glob
import sys

if GRCh38 == "true" :
    ###### Config file and sample sheets #####
    print("Reading configration with GRCh38 version")
    configfile: "config/config_hg38.yaml"
else :
    print("Reading configration with GRCh37 version")
    configfile: "config/config.yaml"


if filetype == "vcf" :
    vcfdata = pd.read_table(config["vcfs"]).set_index("samples", drop=False)
    samplelist = vcfdata.samples.tolist()
else :
    print("filetype cannot be defined")


##### Wildcard constraints #####
if filetype == "vcf" :
    wildcard_constraints:
        sample="|".join(vcfdata.index)
else :
    print("filetype cannot be defined")


if filetype == "vcf" :
    def get_from_vcf(wildcards):
        """Get bam files of given sample-unit."""
        vcfs = vcfdata.loc[(wildcards.sample), "normal"]
        return {"vcf": vcfs}

rule libcrypto_copy:
    output: touch("results/libcrypto.so.1.0.0_copy.done")
    priority: 10000
    log:
        "logs/libcryotp_copy.log",
    shell: "shell/libcrypto_copy.sh 2> {log}"
