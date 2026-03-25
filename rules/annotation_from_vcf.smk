rule vcf2maf:
    input:
        unpack(get_from_vcf)
    output:
        outmaf="results/08vcf2maf/{sample}/{sample}_somatic_filtered_tommo_sed_vep.maf",
        outvcf="results/08vcf2maf/{sample}/{sample}_somatic_filtered_tommo_sed.vcf",
    params:
        ref = config["reference"],
        cachedir = config["cachedir"],
        tommo = config["tommo"],
        clinvar = config["clinvar"],
        spliceai_snv = config["spliceai_snv"],
        spliceai_indel = config["spliceai_indel"],
    conda:
        "../envs/vcf2maf.yaml"
    threads: 16
    log:
        "logs/08vcf2maf/{sample}.log",
    shell:
        "shell/annotation.sh {input.vcf} {output.outmaf} {GRCh38} {params.ref} {params.cachedir} {wildcards.sample} {params.tommo} {output.outvcf} {params.clinvar} {params.spliceai_snv} {params.spliceai_indel} 2> {log}"

rule oncokb_annotator:
    input:
        inputmaf="results/08vcf2maf/{sample}/{sample}_somatic_filtered_tommo_sed_vep.maf",
    output:
        outputmaf="results/09oncokb-annotator/{sample}/{sample}_vep_oncoanno_raw.maf",
    params:
        OncoDir = config["OncoDir"],
        token = config["token"],
    conda:
        "../envs/annotation.yaml"
    threads: 16
    log:
        "logs/09oncokb-annotator/{sample}.log",
    shell:
        "python {params.OncoDir}/MafAnnotator.py -i {input.inputmaf} -o {output.outputmaf} -b {params.token} 2> {log}"

rule intervar:
    input:
        inputvcf="results/08vcf2maf/{sample}/{sample}_somatic_filtered_tommo_sed.vcf",
    output:
        output="results/10intervar/{sample}/{sample}.hg{version}_multianno.txt.intervar",
        output2="results/10intervar/{sample}/{sample}.hg{version}_multianno.txt",
    params:
        intervar = config["intervar"],
        annovar = config["annovar"],
    wildcard_constraints:
        version="\d+"
    conda:
        "../envs/annotation.yaml"
    threads: 16
    log:
        "logs/10intervar/{sample}_hg{version}.log",
    shell:
        "shell/intervar.sh {GRCh38} {params.intervar} {input.inputvcf} {output.output} {params.annovar} {wildcards.sample} {output.output2} 2> {log}"
