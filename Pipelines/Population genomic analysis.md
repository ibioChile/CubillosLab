# Population genomic analysis

This pipeline was used to study the population genomics, phylogeny and structure of Saccharomyces eubayanus, using short-read DNA samples suquenced with Illumina technology.


## Create conda environment to run the programs needed

```
conda create -n snps
conda activate snps
conda install -c bioconda samtools=1.9 freebayes vcftools bowtie2 bcftools bedops bedtools fastqc trimmomatic gatk4
```

## Mapping samples to the genome of the S. eubayanus, CBS12357 strain.

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

## SNP-calling: GATK

1. Index genome file to run [GATK](https://gatk.broadinstitute.org/hc/en-us) for genomic variants detection.

```
gatk CreateSequenceDictionary -R databases/CBS12357_polished_20170509.fasta

samtools faidx databases/CBS12357_polished_20170509.fasta
```

2. Run [GATK]((https://gatk.broadinstitute.org/hc/en-us). Create this script and save it as ```Gatk_CBS12357_scr.sh```.

```
mkdir vcf_gatk_CBS12357

base=${1##*/}; gatk HaplotypeCaller -R databases/CBS12357_polished_20170509.fasta -I $1 --sample-ploidy 2 -O vcf_gatk_CBS12357/${base%.*}.gatk.vcf -ERC GVCF --annotation AlleleFraction --native-pair-hmm-threads 10 --min-base-quality-score 20
```

Now, run (20 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P20 sh Gatk_CBS12357_scr.sh```

The file *VSCBS12357.list* is a file contining a list of *.fixmate.sorted.dedum.bam, separated by a new line. 

3. Join vcf files

```
ls vcf_gatk_CBS12357/*.vcf > gatk_vcf.list

gatk CombineGVCFs -R databases/CBS12357_polished_20170509.fasta -V gatk_vcf.list -O all_samples_combined.gatk.vcf

gatk GenotypeGVCFs -R databases/CBS12357_polished_20170509.fasta --variant all_samples_combined.gatk.vcf --output all_samples_variants.gatk.vcf
```

## Phylogenetic Analysis:

The details of this pipeline can be found [here](https://github.com/carvillarroel/Genetics/blob/master/Genetic%20Analyses%20from%20WGS%20-%20Phylogenetic%20tree%20and%20structure.md)

In summary:

```
vcftools --max-missing 1 --max-alleles 2 --vcf all_samples_variants.gatk.vcf --recode --recode-INFO-all --out GATK_CBS12357_filt_merged
python vcf2phylip -i GATK_CBS12357_filt_merged.recode.vcf
iqtree -s GATK_CBS12357_filt_merged.recode.min4.phy -st DNA -o CL1105.1 -m GTR+ASC -nt 10
iqtree -s GATK_CBS12357_filt_merged.recode.min4.phy.varsites.phy -st DNA -o CL1105.1  -m GTR+ASC -nt 10 -bb 1000 -redo
```

## Structure Analysis:

The details of this pipeline can be found [here](https://github.com/carvillarroel/Genetics/blob/master/Genetic%20Analyses%20from%20WGS%20-%20Phylogenetic%20tree%20and%20structure.md)

In summary:

```
grep -Ev '^CBS12357_mtDNA_polished' all_samples_variants.gatk.vcf > all_samples_variants.gatk.filt.vcf
vcftools --remove remove.txt --vcf all_samples_variants.gatk.filt.vcf --recode --recode-INFO-all --non-ref-ac-any 1 --out onlyeub
./plink_pruning_prep.sh onlyeub.recode.vcf
plink --vcf onlyeub.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.2 --out onlyeub_ldfilter --threads 10

sed 's/polished_/polished\t/' onlyeub_ldfilter.prune.in > onlyeub_ldfilter.prune.in.vcftools
vcftools --vcf onlyeub.recode.vcf --positions onlyeub_ldfilter.prune.in.vcftools --recode-INFO-all --recode --out onlyeub_ldfilter

vcftools --vcf onlyeub_ldfilter.recode.vcf --thin 1000 --recode --recode-INFO-all --out onlyeub_ldfilter_thinned

mkdir onlyeub_structure
populations -t 10 -V onlyeub_ldfilter_thinned.recode.vcf -O onlyeub_structure --structure

cd onlyeub_structure
tail -n +2 onlyeub_ldfilter_thinned.recode.p.structure > onlyeub_structure.txt
cd ..

for m in {2..8}
do
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep1 &$
sleep 30
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep2 &$
sleep 30
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep3 &$
sleep 30
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep4 &$
sleep 30
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep5 &$
sleep 30
sbatch --wrap="structure -m mainparams -p extraparams -K $m -L 10354 -N 285 -i onlyeub_structure/onlyeub_structure.txt -o onlyeub_structure_res_K.${m}_rep6 &$
sleep 30
done
```

2. Run vcfallelic to compare among programs. Save the following script as ``` vcfallelicprim_CBS12357.sh```.

```base=${1##*/}; vcfallelicprimitives vcf_freebayes_CBS12357_filt/${base%_*}.fb.filt.vcf > vcf_freebayes_CBS12357_filt/${base%_*}.fb.filt.all.vcf;  vcfallelicrimitives vcf_gatk_CBS12357_filt/${base%_*}.gatk.filt.vcf > vcf_gatk_CBS12357_filt/${base%_*}.gatk.filt.all.vcf```

Now, run (50 parallel processes):

```cat VSCBS12357.list | xargs -n1 -P50 sh vcfallelicprim_CBS12357.sh```

3. Compare SNPs in all strains using the pyth
