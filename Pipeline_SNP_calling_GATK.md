# Pipeline for SNP calling using GATK

PREPROCESS READS

>batch processes samples in FASTQC to assess overall quality, per read quality, and presence of adaptor sequence
fastqc --threads 4 [FILE1] [FILE2] [FILE n]

>get read quality sumary report
qualimap .

>Use fastp for trimming low quality bases (3'), discard low quality reads (too many Ns, length < 36 bp), and trimming adaptors en SRA sequence files

fastp --in1 9161196_S318_R1_001_paired.fastq.gz --in2 9161196_S318_R2_001_paired.fastq.gz --out1 705.1_Talca_1.fastq.gz --out2 705.1_Talca_2.fastq.gz -A -3 -l 36 -w 4 -j 705.1_Talca.json -h 705.1_Talca.html

MAP TO REFERENCE

>map cleaned reads to reference using BWA-mem (include sample information as read group @RG, mark secondary alignments -M)
bwa index CBS12357_polished_20170509.fasta
bwa mem -t 4 -M -R '@RG\tID:701.1_Talca\tSM:701.1_Talca' -o 701.1_Talca_al.sam CBS12357_polished_20170509.fasta   701.1_Talca_1.fastq.gz  701.1_Talca_2.fastq.gz

>sort  alignments (use bam output)
samtools sort --threads 4 -o 701.1_Talca_sorted.sam 701.1_Talca_al.sam &&

>Mark duplicates 
picard MarkDuplicates I= 701.1_Talca_sorted.sam O= 701.1_Talca_sorted_DP.sam M= 701.1_Talca__DP.txt &&

>Index alignments
samtools index 701.1_Talca_sorted_DP.sam

>get mapping quality report
qualimap bamqc -c -bam 1007.1_TDP_sorted_DP.bam -outformat HTML -outdir /mapping_reports/1007.1_TDP_sorted_DP_qualimap_report

>get summary report (in mapping_reports folder)
multiqc .


CALL VARIANTS
>build dictionary file for reference fasta
gatk CreateSequenceDictionary -R CBS12357_polished_20170509.fasta

>Intervals each chromosome Eubayanus ref genome (nanopore)
--intervals CBS12357_Chr01_polished:1-210599 --intervals CBS12357_Chr02_polished:1-1276473 --intervals CBS12357_Chr03_polished:1-308584 --intervals CBS12357_Chr04_polished:1-985108 --intervals CBS12357_Chr05_polished:1-586080 --intervals CBS12357_Chr06_polished:1-263164 --intervals CBS12357_Chr07_polished:1-1057760 --intervals CBS12357_Chr08_polished:1-834776 --intervals CBS12357_Chr09_polished:1-411897 --intervals CBS12357_Chr10_polished:1-752744 --intervals CBS12357_Chr11_polished:1-643305 --intervals CBS12357_Chr12_polished:1-1150293 --intervals CBS12357_Chr13_polished:1-958629 --intervals CBS12357_Chr14_polished:1-767356 --intervals CBS12357_Chr15_polished:1-742518 --intervals CBS12357_Chr16_polished:1-914559 --intervals CBS12357_mtDNA_polished:1-86031 -


>Call variants per sample and chromosome

gatk HaplotypeCaller --intervals CBS12357_mtDNA_polished:1-86031 --emit-ref-confidence GVCF -R CBS12357_polished_20170509.fasta    -I 1001.1_TDP_sorted_DP.bam    -O 1001.1_TDP_sorted_mtDNA.g.vcf 


>build variant db for gatk

gatk GenomicsDBImport --genomicsdb-workspace-path all_variants_chr1/ --intervals CBS12357_Chr01_polished:1-210599   --arguments_file samples1.txt&


>call genotypes from all samples

gatk GenotypeGVCFs -R CBS12357_polished_20170509.fasta -V gendb://all_variants_chr1/  -G StandardAnnotation -O variants_chr1.vcf 2>chr1_log.txt&

>merge variant files generated for each chromosome 

gatk MergeVcfs -I variant_all_files.list -O all_raw_variants2.vcf

>quality "hard filters" for SNPs and INDELs

gatk SelectVariants -R CBS12357_polished_20170509.fasta -V all_raw_variants2.vcf --select-type-to-include SNP -O raw_snps.vcf

gatk VariantFiltration -R CBS12357_polished_20170509.fasta -V raw_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Se_snp_filter" -O filtered_snps.vcf

gatk SelectVariants -R CBS12357_polished_20170509.fasta -V all_raw_variants2.vcf --select-type-to-include INDEL -O raw_INDELs.vcf

gatk VariantFiltration -R CBS12357_polished_20170509.fasta -V raw_INDELs.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "Se_indel_filter" -O filtered_INDELs.vcf 




