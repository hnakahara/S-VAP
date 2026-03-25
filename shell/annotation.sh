#!/bin/sh

input=$1
output=$2
GRCh38=$3
ref=$4
cachedir=$5
tumor=$6
tommo=$7
outvcf=$8
veppath=$(which vep)
tmp1=${input%.vcf.gz}_tommo.vcf
tmp2=${input%.vcf.gz}_tommo_clin.vcf
clinvar=$9
spliceai_snv=$10
spliceai_indel=$11

if $GRCh38 ; then
  echo "GRCh38: $GRCh38"
  echo "build version GRCh38"
  echo "Start adding ToMMo"
  SnpSift '-Xmx24g' annotate $tommo $input > $tmp1
  SnpSift '-Xmx24g' annotate $clinvar $tmp1 > $tmp2
  #echo "replacement AF to AF_ToMMo"
  sed 's/AF=/AF=/g' $tmp2 > $outvcf
  echo "Start vcf2maf"
  vcf2maf.pl --input-vcf $outvcf --output-maf $output --ref-fasta $ref --vep-data $cachedir --vep-path ${veppath%/*} --vep-forks 24 --species homo_sapiens --ncbi-build GRCh38 --cache-version 112 --tumor-id $tumor --vep-custom $tommo,ToMMo,vcf,exact,,AF --vep-plugins SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel},cutoff=0.5 --retain-ann ToMMo_AF,SpliceAI_cutoff,SpliceAI_pred --retain-info AC_XX,AC_XY,AN_XX,AN_XY,NAME,Report,MAF,CLNSIG,CLNSIGCONF,CLNREVSTAT --max-subpop-af 0.001
  rm -rf $tmp1
  rm -rf $tmp2
else
  echo "GRCh38: $GRCh38"
  echo "build version GRCh37"
  echo "Start adding ToMMo"
  SnpSift '-Xmx24g' annotate $tommo $input > $tmp1
  SnpSift '-Xmx24g' annotate $clinvar $tmp1 > $tmp2
  #echo "replacement AF to AF_ToMMo"
  sed 's/AF=/AF=/g' $tmp2 > $outvcf
  echo "Start vcf2maf"
  vcf2maf.pl --input-vcf $outvcf --output-maf $output --ref-fasta $ref --vep-data $cachedir --vep-path ${veppath%/*} --vep-forks 24 --species homo_sapiens --ncbi-build GRCh37 --cache-version 112 --tumor-id $tumor --vep-custom $tommo,ToMMo,vcf,exact,,AF --vep-plugins SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel},cutoff=0.5 --retain-ann ToMMo_AF,SpliceAI_cutoff,SpliceAI_pred --retain-info AC_XX,AC_XY,AN_XX,AN_XY,NAME,Report,MAF,CLNSIG,CLNSIGCONF,CLNREVSTAT --max-subpop-af 0.001
  rm -rf $tmp1
  rm -rf $tmp2
fi
