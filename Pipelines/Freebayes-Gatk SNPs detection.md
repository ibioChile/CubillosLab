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

2. Select variants unique to each strain.

```
bgzip -c 712B_N.crassa.vcf > 712B_N.crassa.vcf.gz
tabix -f -p vcf 712B_N.crassa.vcf.gz

bgzip -c 712A_N.crassa.vcf > 712A_N.crassa.vcf.gz
tabix -f -p vcf 712A_N.crassa.vcf.gz

perl /bin/vcftools/src/perl/vcf-isec -c 712A_N.crassa.vcf.gz 712B_N.crassa.vcf.gz > 712A_N.unique.vcf

perl /bin/vcftools/src/perl/vcf-isec -c 712B_N.crassa.vcf.gz 712A_N.crassa.vcf.gz > 712B_N.unique.vcf
```

3. Filter variants by coverage (> 10 mapping reads) and quality (> 30).

```
bcftools view -i 'FORMAT/DP>10' 712B_N.unique.vcf > 712B_N.unique_filtDP.vcf 
bcftools view -i 'FORMAT/DP>10' 712A_N.unique.vcf > 712A_N.unique_filtDP.vcf 

vcftools --vcf 712A_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712A_N.unique_filtered
vcftools --vcf 712B_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712B_N.unique_filtered
```

## GATK

1. Index genome file to run [GATK](https://gatk.broadinstitute.org/hc/en-us) for genomic variants detection.

```
gatk CreateSequenceDictionary -R databases/CBS12357_polished_20170509.fasta

samtools faidx databases/CBS12357_polished_20170509.fasta
```

2. Run [GATK]((https://gatk.broadinstitute.org/hc/en-us). Create this script and save it as ```Gatk_CBS12357_scr.sh```.

```
gatk HaplotypeCaller -R GCF_000182925.2_NC12_genomic.fna -I 712A_N.crassa.fixmate.sorted.dedup.RG.bam --sample-ploidy 1 -O 712A_N.crassa.vcf

gatk HaplotypeCaller -R GCF_000182925.2_NC12_genomic.fna -I 712B_N.crassa.fixmate.sorted.dedup.RG.bam --sample-ploidy 1 -O 712B_N.crassa.vcf
```

```base=${1##*/}; gatk HaplotypeCaller -R databases/CBS12357_polished_20170509.fasta -I $1 --sample-ploidy 2 -O vcf_gatk_CBS12357/${base%.*}.gatk.vcf --annotation AlleleFraction --native-pair-hmm-threads 10 --min-base-quality-score 20```

Now, run (20 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P20 sh Gatk_CBS12357_scr.sh```


4. Select variants unique to each strain.

```
bgzip -c 712B_N.crassa.vcf > 712B_N.crassa.vcf.gz
tabix -f -p vcf 712B_N.crassa.vcf.gz

bgzip -c 712A_N.crassa.vcf > 712A_N.crassa.vcf.gz
tabix -f -p vcf 712A_N.crassa.vcf.gz

perl /bin/vcftools/src/perl/vcf-isec -c 712A_N.crassa.vcf.gz 712B_N.crassa.vcf.gz > 712A_N.unique.vcf

perl /bin/vcftools/src/perl/vcf-isec -c 712B_N.crassa.vcf.gz  712A_N.crassa.vcf.gz > 712B_N.unique.vcf
```

5. Filter variants by coverage (> 10 mapping reads) and quality (> 30).

```
bcftools view -i 'FORMAT/DP>10' 712B_N.unique.vcf > 712B_N.unique_filtDP.vcf 
bcftools view -i 'FORMAT/DP>10' 712A_N.unique.vcf > 712A_N.unique_filtDP.vcf 

vcftools --vcf 712A_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712A_N.unique_filtered
vcftools --vcf 712B_N.unique_filtDP.vcf  --recode --recode-INFO-all --minQ 30 --out 712B_N.unique_filtered
```
