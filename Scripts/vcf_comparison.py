import vcf
import glob
import os
import pandas as pd
from multiprocessing.pool import Pool

arr = glob.glob("*216_filt/*.all.vcf")

snps = []

def vcf_fun(vcf_file):
        print("working on file "+ vcf_file)
        sample = os.path.split(vcf_file)[1].split("_")[0]
        reference = os.path.split(vcf_file)[1].split("_")[1]
        program = os.path.split(vcf_file)[1].split("DP.")[1].split(".")[0]
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        snps=[]
        for record in vcf_reader:
                snps.append([sample, program, reference, record.CHROM, str(record.POS), str(record.REF), str(record.ALT[0])])
        return snps

if __name__ == '__main__':
    p = Pool(30)
    snps = p.map(vcf_fun, [vcf_file for vcf_file in arr])
    snps = [ent for sublist in snps for ent in sublist]

print("creating table")
df = pd.DataFrame(snps, columns=['Sample','Program','Reference strain','CHROM','POS','REF','ALT'])
print("creating pivot table")
df_pivot = df.pivot_table(index=['CHROM','POS','REF','ALT'], columns=['Sample','Program','Reference strain'], aggfunc=len).fillna(0)
df_pivot.to_csv(r'comparison_SNPs_216_programs.csv', index=True)
