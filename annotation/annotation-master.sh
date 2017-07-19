#!/usr/bin/bash

if [ $# -lt 2 ]
then
    echo "\nUsage:\n\tbash $0 <variants[vcf]> <output-dir>\n"
    exit 0
fi

vcf=$1
wd=$2
transcriptList=/Users/mmmskje2/Desktop/intronics/by-chr/RD_v3_transcripts.txt
cd $wd
mkdir $wd/tmp
    ##  filter vcf

chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

    ##  vcf --> vep --> filters [check numbers]
mkdir $wd/tmp/vep
echo "new directory created: $wd/tmp/vep"
echo "using docker to annotate with VEP"
docker run -v /Users/mmmskje2/.vep/:/home/vep/.vep -v /Users/mmmskje2/.vep/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz:/home/vep/.vep/in.fa -v $wd/tmp/vep:/home/vep/src/ensembl-vep/output/ -v $vcf:/home/vep/src/ensembl-vep/input/in.txt --rm willmclaren/ensembl-vep ./vep -i input/in.txt -o output/out-refseq.vep --offline --force_overwrite --assembly GRCh37 --everything --refseq --hgvs --fasta /home/vep/.vep/in.fa
echo "selecting transcripts"
perl ~/Desktop/intronics/by-chr/filter-vep-file.pl --vep $wd/tmp/vep/out-refseq.vep --outputdir $wd/tmp/vep --transcripts $transcriptList



    ##  vcf --> gnomad [check numbers]
mkdir $wd/tmp/gnomad
echo "new directory created: $wd/tmp/gnomad"
echo "annotating with gnomad frequencies"
perl ~/Documents/git/gnomad/get_gnomad_freqs_withIndels.pl --variants $vcf --outputdir $wd/tmp/gnomad --type vcf

    ##  vcf --> inhouse [check numbers]
### not suitable for macbook

mkdir $wd/tmp/combine
echo "new directory created: $wd/tmp/combine"
echo "combining vep and gnomad annotations"
perl ~/Documents/git/annotation/combine-vep-gnomad.pl --gnomad $wd/tmp/gnomad/gnomadfreq.txt --vep $wd/tmp/vep/vep-filtered-transcripts.txt --outputdir $wd/tmp/combine

