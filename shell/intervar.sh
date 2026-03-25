#!/bin/bash

GRCh38=$1
intervar=$2
input=$3
output=$4
annovar=$5
sample=$6
dir=$(dirname $output)


if $GRCh38 ; then
  echo "GRCh38: $GRCh38"
  echo "build version GRCh38"
  echo "Start InterVar annotation"
  $intervar/Intervar.py -b hg38 -t $intervar/intervardb -i $input --input_type=VCF -o $dir/$sample --table_annovar=$annovar/table_annovar.pl --convert2annovar=$annovar/convert2annovar.pl --annotate_variation=$annovar/annotate_variation.pl --database_locat=$annovar/humandb
else
  echo "GRCh38: $GRCh38"
  echo "build version GRCh37"
  echo "Start InterVar annotation"
  $intervar/Intervar.py -b hg19 -t $intervar/intervardb -i $input --input_type=VCF -o $dir/$sample --table_annovar=$annovar/table_annovar.pl --convert2annovar=$annovar/convert2annovar.pl --annotate_variation=$annovar/annotate_variation.pl --database_locat=$annovar/humandb
fi
