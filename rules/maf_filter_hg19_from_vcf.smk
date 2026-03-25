rule maf_filter:
    input:
        inputmaf="results/09oncokb-annotator/{sample}/{sample}_vep_oncoanno_raw.maf",
        input_intervar="results/10intervar/{sample}/{sample}.hg19_multianno.txt.intervar",
        input_intervar_raw="results/10intervar/{sample}/{sample}.hg19_multianno.txt"
    output:
        outmaf1="results/11final-output/{sample}/{sample}_variant_class.maf",
        outmaf2="results/11final-output/{sample}/{sample}_variant_class_vus_prior.maf",
        outmaf3="results/11final-output/{sample}/{sample}_merged_raw.maf",
    params:
        #genefile = config["genefile"],
    #wildcard_constraints:
    #    version="\d+"
    conda:
        "../envs/annotation.yaml"
    threads: 16
    log:
        "logs/11final-output/maffilter/{sample}.log",
    shell:
        "Rscript Rscript/maf_filter.r {input.inputmaf} {output.outmaf1} {output.outmaf2} {output.outmaf3} {input.input_intervar} {input.input_intervar_raw} 2> {log}"


rule maf_merge:
    input:
        inputmaf1=expand("results/11final-output/{sample}/{sample}_variant_class.maf", sample=samplelist),
        inputmaf2=expand("results/11final-output/{sample}/{sample}_variant_class_vus_prior.maf", sample=samplelist),
        inputmaf3=expand("results/11final-output/{sample}/{sample}_merged_raw.maf", sample=samplelist),
    output:
        outmaf1="results/11final-output/merge_variant_class.maf",
        outmaf2="results/11final-output/merge_variant_class_vus_prior.maf",
        outmaf3="results/11final-output/merge_raw.maf",
    params:
        dir="results/11final-output"
    conda:
        "../envs/annotation.yaml"
    threads: 24
    log:
        "logs/11final-output/maffilter/maf_merge.log",
    shell:
        "Rscript Rscript/maf_merge.r {params.dir} {output.outmaf1} {output.outmaf2} {output.outmaf3} 2>> {log}"
