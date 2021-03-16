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

Optimal K values can be obtained using [structure-selector](http://lmme.qdio.ac.cn/ StructureSelector/) (Li & Liu, 2018) according to the Evanno method (Evanno et al. 2005). The resulting plots were obtained using CLUMPAK (Kopelman et al. 2015) and visualized using [Structure Plot V2.0](http://omicsspeaks.com/strplot2/) (Ramasamy et al. 2014). 

## Admixture Analysis:

Install Admixture:

```conda install -c bioconda admixture```

Script:

```
plink --vcf onlyeub.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 5 0.2 --maf 0.05 --out onlyeub_ldfilter --make-bed --threads 10
awk '{$1=0;print $0}' onlyeub_ldfilter.bim > onlyeub_ldfilter.bim.tmp
mv  onlyeub_ldfilter.bim.tmp  onlyeub_ldfilter.bim

for k in {2..10}
do
admixture -j20 --cv onlyeub_ldfilter.bed $k > onlyeub_ldfilter.$k.log
done
```

Use [structure-selector](http://lmme.qdio.ac.cn/ StructureSelector/) (Li & Liu, 2018) to select K value.

## FineStructure Analysis:

The details of this pipeline can be found [here](https://github.com/carvillarroel/Genetics/blob/master/Genetic%20Analyses%20from%20WGS%20-%20fineStructure%20and%20Globetrotter.md)

In summary:

```
SnpSift split onlyeub.recode.vcf

plink --vcf onlyeub.recode.CBS12357_Chr01_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr1
plink --vcf onlyeub.recode.CBS12357_Chr02_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr2
plink --vcf onlyeub.recode.CBS12357_Chr03_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr3
plink --vcf onlyeub.recode.CBS12357_Chr04_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr4
plink --vcf onlyeub.recode.CBS12357_Chr05_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr5
plink --vcf onlyeub.recode.CBS12357_Chr06_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr6
plink --vcf onlyeub.recode.CBS12357_Chr07_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr7
plink --vcf onlyeub.recode.CBS12357_Chr08_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr8
plink --vcf onlyeub.recode.CBS12357_Chr09_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr9
plink --vcf onlyeub.recode.CBS12357_Chr10_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr10
plink --vcf onlyeub.recode.CBS12357_Chr11_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr11
plink --vcf onlyeub.recode.CBS12357_Chr12_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr12
plink --vcf onlyeub.recode.CBS12357_Chr13_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr13
plink --vcf onlyeub.recode.CBS12357_Chr14_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr14
plink --vcf onlyeub.recode.CBS12357_Chr15_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr15
plink --vcf onlyeub.recode.CBS12357_Chr16_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr16

for i in {1..16}
do
bash phasing_pipeline/run.sh chr${i}.ped chr${i}.map chr${i}
done

for i in {1..16}
do
perl fs_4.1.1/plink2chromopainter.pl -p=chr${i}.phased.ped -m=chr${i}.phased.map -o=chr${i}.chromopainter -f
done

for i in {1..16}
do
perl fs_4.1.1/makeuniformrecfile.pl chr${i}.chromopainter chr${i}_rec
done

for i in {1..2};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {3..4};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {5..6};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {7..8};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {9..10};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {11..12};do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {13..14; do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
for i in {15..16}; do fs_4.1.1/fs_linux_glibc2.3 chromopainter -g chr${i}.chromopainter -r chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done

fs_4.1.1/fs_linux_glibc2.3 chromocombine -d datos
fs_4.1.1/fs_linux_glibc2.3 finestructure  -x 100000 -y 100000 -z 1000 -X -Y output.chunkcounts.out out.105eubs.mcmc.xml
fs_4.1.1/fs_linux_glibc2.3 finestructure -x 100000 -k 2 -m T -t 1000000 -X -Y output.chunkcounts.out out.105eubs.mcmc.xml structure_tree.out
fs_4.1.1/fs_linux_glibc2.3 finestructure -X -Y -e meancoincidence output.chunkcounts.out out.105eubs.mcmc.xml structure_meancoincidence.csv
fs_4.1.1/fs_linux_glibc2.3 finestructure -X -Y -e X2  output.chunkcounts.out out.105eubs.mcmc.xml structure_meanstate.csv
```
