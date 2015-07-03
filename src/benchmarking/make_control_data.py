""" generates test data for mosaic de novo calling
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import subprocess
import os

import pysam

BENCHMARKING_DIR = "/lustre/scratch113/projects/ddd/users/jm33/mosaic_benchmarking/"
SAMTOOLS = "/software/vertres/bin-external/samtools-0.1.18"
BEDTOOLS = "/software/vertres/bin-external/bedtools"
FIX_READ_GROUP = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/src/benchmarking/fix_read_group.py"

# # copy all of Art's files
# rsync -vaP /lustre/scratch113/projects/mu/users/aw15/samtoolsMod/* /lustre/scratch113/projects/ddd/users/jm33/mosaic_benchmarking
# # prepare a temp dir
# mkdir /lustre/scratch113/projects/ddd/users/jm33/mosaic_benchmarking/temp
# cp /lustre/scratch113/projects/ddd/users/jm33/mosaic_benchmarking/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam.21.bam /lustre/scratch113/projects/ddd/users/jm33/mosaic_benchmarking/temp/temp.21.bam

def get_median_coverage(bam_path, chrom):
    """ get the median coverage for a BAM on a specific chromosome
    
    Args:
        bam_path: path to BAM file
        chrom: chromosome that we wish to obtain the coverage on
    
    Returns:
        median coverage as an integer
    """
    
    # get the median coverage across one chrom of the BAM
    samtools = subprocess.Popen([SAMTOOLS, "view", "-b", \
        os.path.join(BENCHMARKING_DIR, bam_path), chrom], stdout=subprocess.PIPE)
    coverage = subprocess.check_output([BEDTOOLS, "genomecov", "-ibam", \
        "stdin"], stdin=samtools.stdout)
    
    coverage = coverage.decode("utf-8").split("\n")
    
    adjusted_chrom_length = None
    zero_depth_bases = None
    cumulative_prob = 0
    median_coverage = 0
    for line in coverage:
        line = line.strip().split()
        current_chrom = line[0]
        depth = int(line[1])
        n_bases = int(line[2])
        chrom_length = int(line[3])
        
        if current_chrom != chrom:
            continue
        
        # if we are looking at the bases with zero depth, use this to adjust the
        # chrom length to exclude the unsequenceable regions of the chrom. This
        # doesn't adjust the median coverage by much, just bumps up by 2-4X.
        if zero_depth_bases == None and depth == 0:
            zero_depth_bases = n_bases
        
        # when we calculate the cumulative proportion, we need to account for
        # the unsequenceable section of the chromosome, which is given by the
        # bases with zero depth
        if adjusted_chrom_length == None:
            adjusted_chrom_length = chrom_length - zero_depth_bases
        
        if depth > 0:
            cumulative_prob += n_bases/adjusted_chrom_length
        
        # adjust the median coverage, unless we are past the median
        if cumulative_prob < 0.5:
            median_coverage = depth
        else:
            break
    
    return median_coverage

def subsample_bam(input_bam, alternate_bam, merged_path, proportion):
    """ subsample a BAM to get a reduced proportion of reads
    
    Args:
        input_bam: path to BAM to subsample
        alternate_bam: path to BAM to merge with the subsampled BAM
        merged_path: path to BAM to output merged BAM to.
        proportion: proportion of the input BAM reads to retain
    
    Returns:
        subsampled BAM as a subprocess PIPE.
    """
    
    assert proportion > 0
    assert proportion < 1
    assert os.path.exists(input_bam)
    assert os.path.exists(alternate_bam)
    
    # samtools can subsample, e.g. to create a BAM with 10% of the initial reads
    probability = "1." + "{}".format(proportion)[2:]
    subsample = subprocess.Popen([SAMTOOLS, "view", "-b", "-s", probability, \
        input_bam], stdout=subprocess.PIPE)
    
    # merge the bams with: samtools merge MERGED_PATH BAM1_PATH BAM2_PATH
    # where MERGED_PATH can be "-" for standard out, and one BAM_PATH can be
    # /dev/stdin for stdin (merge cannot use the normal stdin pipe).
    merge = subprocess.Popen([SAMTOOLS, "merge", "-", alternate_bam, \
        "/dev/stdin"], stdin=subsampled.stdout, stdout=subprocess.PIPE)
    
    # and fix the read group information in the merged BAM
    subprocess.call(["python3", FIX_READ_GROUP, \
        "--subsampled", bam1_path, "--alternate", bam2_path], \
        stdin=merged.stdout, stdout=open(merged_path, "w"))

def main():
    
    bam1_path = os.path.join(BENCHMARKING_DIR, "CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam.21.bam")
    bam2_path = os.path.join(BENCHMARKING_DIR, "CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam.21.bam")
    merged_bam = os.path.join(BENCHMARKING_DIR, "TEMP.bam")
    
    # bam1_coverage = get_median_coverage(bam1_path, "21")
    # bam2_coverage = get_median_coverage(bam2_path, "21")
    
    # proportion = 0.5 /((0.5 * bam_1_coverage) / bam2_coverage)
    proportion = 0.8
        
    subsample_bam(bam1_path, bam2_path, merged_bam, proportion)

    

if __name__ == '__main__':
    main()




# # then run modified and standard samtools to call de novos
# SAMTOOLS mpileup -p0.5 -uf ref.fa sample.bam | bcftools view -vcg - > var.raw.vcf
