

samtools mpileup -gDf reference.fa father.bam > father.bcf
samtools mpileup -gDf reference.fa mother.bam > mother.bcf
 
# Create a BCF for the child using the modified Samtools. In this case, I assume that you use /nfs/team29/aw15/samtoolsMod/samtools, and a probability of 25% that a read supports de novo mutations.
 
/nfs/team29/aw15/samtoolsMod/samtools mpileup -p0.25 -gDf reference.fa child.bam > child.bcf
 
# Convert all three BCFs to VCFs:
 
bcftools view father.bcf > father.vcf
bcftools view mother.bcf > mother.vcf
bcftools view child.bcf > child.vcf
 
# Prepare VCFs for merging by zipping and tabix indexing:
 
bgzip father.vcf
bgzip mother.vcf
bgzip child.vcf

tabix father.vcf.gz
tabix mother.vcf.gz
tabix child.vcf.gz
 
# Do the actual merge using Vcftools:
 
/nfs/team29/aw15/software/vcftools_0.1.9/bin/vcf-merge father.vcf.gz mother.vcf.gz child.vcf.gz > merged.vcf.temp
 
# Fix the VCF PL fields. This step is necessary because when you merged the three VCFs, in some cases there may not be enough information to determine genotype likelihoods for genotypes involving an alternative allele in all three samples. In such a case Vcftools' merge function substitutes the genotype likelihood in the resulting VCF's PL field with a dot ("."). This is not something DNG can work with, and therefore the dots have to be substituted with extremely low phred-scaled genotype likelihoods such as 255. This is what the plFix.pl script does. Please note that this script is quite crude and you should probably understand what it does before using it for something important
 
perl /nfs/team29/aw15/scripts/plFix.pl merged.vcf.temp > merged.vcf
 
# Convert the VCF to BCF. seqDic.txt is a file that simply contains the name of each chromosome in your VCF in a separate line, e.g. 1, 2, 3, …, 22, X, Y. Or: chr1, chr2, chr3, …, chr22, chrX, chrY.
 
bcftools view -D seqDic.txt -Sb merged.vcf > merged.bcf
 
# Use DNG to call de novo mutations. In this particular case, I use v0.5. As always with DNG, the appropriate PED file (family.ped) has to be in the working directory.
 
/nfs/users/nfs_s/sa9/scripts/denovogear/denovogear-0.5/build/src/denovogear dnm auto --ped family.ped --bcf merged.bcf > merged.dnm
 


