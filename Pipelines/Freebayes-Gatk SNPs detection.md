# Freebayes-GATK SNPs detection

This pipeline uses short sequencing reads to detect genomic variants (insertions, deletions and SNPs) in genome. Here, we compare the output of two known variant callers, freebayes and GATK. This specific example evaluates variants in strains of Saccharomyces eubayanus, when compared to the strain CBS12357. 

## Create conda environment to run the programs needed

```
conda create -n snps
conda activate snps
conda install -c bioconda samtools=1.9 freebayes vcftools bowtie2 bcftools bedops bedtools fastqc trimmomatic gatk4
```

## Fix mate in bam files, remove duplicates and sort.

```
for file in *.bam;
do
base=${file##*/}

if [ -s ${base%.*}.fixmate.sorted.dedup.bam.bai ]
then
    	echo "$base already processed"
else
    	  samtools sort -n -o ${base%.*}.sorted.bam $file
        samtools fixmate -m ${base%.*}.sorted.bam ${base%.*}.fixmate.bam
        samtools sort ${base%.*}.fixmate.bam -o ${base%.*}.fixmate.sorted.bam
        samtools markdup -r ${base%.*}.fixmate.sorted.bam ${base%.*}.fixmate.sorted.dedup.bam
        samtools index ${base%.*}.fixmate.sorted.dedup.bam
        rm ${base%.*}.sorted.bam ${base%.*}.fixmate.bam ${base%.*}.fixmate.sorted.bam
fi
done
```

## Freebayes

1. Run [freebayes](https://github.com/ekg/freebayes) for genomic variants detection. Create this script and save it as ```Freebayes_CBS12357_scr.sh```.

```base=${1##*/}; freebayes -f databases/CBS12357_polished_20170509.fasta -p 2 $1 > vcf_freebayes_CBS12357/${base%.*}.fb.vcf```

Now, run (30 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P30 sh Freebayes_CBS12357_scr.sh```

Where ```VSCBS12357.list``` lists all fixed bam files paths.


## GATK

1. Index genome file to run [GATK](https://gatk.broadinstitute.org/hc/en-us) for genomic variants detection.

```
gatk CreateSequenceDictionary -R databases/CBS12357_polished_20170509.fasta

samtools faidx databases/CBS12357_polished_20170509.fasta
```

2. Run [GATK]((https://gatk.broadinstitute.org/hc/en-us). Create this script and save it as ```Gatk_CBS12357_scr.sh```.

```base=${1##*/}; gatk HaplotypeCaller -R databases/CBS12357_polished_20170509.fasta -I $1 --sample-ploidy 2 -O vcf_gatk_CBS12357/${base%.*}.gatk.vcf --annotation AlleleFraction --native-pair-hmm-threads 10 --min-base-quality-score 20```

Now, run (20 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P20 sh Gatk_CBS12357_scr.sh```


## Variants filter and comparison

1. Filter variants by coverage > 4 mapping reads, mapping quality > 30 and phred score > 20. Save the following script as ```bcftools_CBS12357_scr.sh```.

```
base=${1##*/}; bcftools view -i 'FORMAT/AO>4 & MQM>30 & FORMAT/QA>20' vcf_freebayes_CBS12357/${base%.*}.fb.vcf  > vcf_freebayes_CBS12357_filt/${base%_*}.fb.filt.vcf; bcftools view -i 'FORMAT/AD>4 & MQ>30' vcf_gatk_CBS12357/${base%.*}.fb.vcf  > vcf_gatk_CBS12357_filt/${base%_*}.fb.filt.vcf
```

Now, run (50 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P50 sh bcftools_CBS12357_scr.sh```

2. Run vcfallelic to compare among programs. Save the following script as ``` vcfallelicprim_CBS12357.sh```.

```base=${1##*/}; vcfallelicprimitives vcf_freebayes_CBS12357_filt/${base%_*}.fb.filt.vcf > vcf_freebayes_CBS12357_filt/${base%_*}.fb.filt.all.vcf;  vcfallelicrimitives vcf_gatk_CBS12357_filt/${base%_*}.gatk.filt.vcf > vcf_gatk_CBS12357_filt/${base%_*}.gatk.filt.all.vcf```

Now, run (50 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P50 sh vcfallelicprim_CBS12357.sh```

