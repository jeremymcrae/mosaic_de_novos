""" a bunch of functions to help with mosaic de novo calling
"""

from __future__ import absolute_import

import os
import pysam

def make_corrected_vcf_header(bam_path, output=None):
    """ makes a header file that fixes the lack of explanatory lines
    
    Args:
        bam_path: path to the bam file that needs a corrected VCF header.
        output: potential file handle for the header file
    """
    
    # here's a VCF header that includes the VCF output types from samtools
    header = ["##fileformat=VCFv4.1\n",
        "##samtoolsVersion=0.1.18 (r982:295)\n",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n",
        "##INFO=<ID=I16,Number=.,Type=Integer,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">\n",
        "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n",
        "##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better).\">\n",
        "##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\">\n",
        "##INFO=<ID=INDEL,Number=1,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n",
        "##INFO=<ID=IS,Number=.,Type=Float,Description=\"Indel score?\">\n",
        "##INFO=<ID=QS,Number=.,Type=Float,Description=\"Auxiliary tag used for calling.\">\n",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"]
    
    bam_dir = os.path.dirname(bam_path)
    if output is None:
        new_path = os.path.join(bam_dir, "fixed_header.txt")
        output = open(new_path, "w")
    
    alt_id = get_sample_id_from_bam(bam_path)
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + alt_id + "\n")
    
    # write the header
    output.writelines(header)
    output.flush()

def make_seq_dic_file(seq_dic):
    """ make sure we have a contig dictionary file available for bcftools
    
    Args:
        seq_dic: path to the sequence dictionary file, or a file handle.
        
    Returns:
        nothing
    """
    
    chroms = list(range(1, 22)) + ["X", "Y"]
    chroms = [str(x) for x in chroms]
    chroms = "\n".join(chroms) + "\n"
    
    if isinstance(seq_dic, str):
        # only create the file if it doesn't already exist
        if os.path.exists(seq_dic):
            return
        
        with open(seq_dic, "w") as output:
            output.writelines(chroms)
    else:
        seq_dic.writelines(chroms)
        seq_dic.flush()

def make_ped_for_trio(child, mother, father, sex, ped):
    """ make a PED file to define the relationship of samples in a BCF file
    
    Args:
        child: path to a bam file for the child.
        mother: path to a bam file for the mother.
        father: path to a bam file for the father.
        ped: path to write the PED file to, or a file handle.
    
    Returns:
        nothing
    """
    
    # get the alternative IDs that will exist in the VCF file
    child_id = get_sample_id_from_bam(child)
    mother_id = get_sample_id_from_bam(mother)
    father_id = get_sample_id_from_bam(father)
    
    # make a ped file
    fam_id = "temp"
    
    # make the lines for the PED file, in PED format.
    child_line = "\t".join([fam_id, child_id, father_id, mother_id, sex, "2"]) + "\n"
    father_line = "\t".join([fam_id, father_id, "0", "0", "1", "1"]) + "\n"
    mother_line = "\t".join([fam_id, mother_id, "0", "0", "2", "1"]) + "\n"
    lines = [child_line, father_line, mother_line]
    
    if isinstance(ped, str):
        # only create the file if it doesn't already exist
        if os.path.exists(ped):
            return
        
        with open(ped, "w") as output:
            output.writelines(lines)
    else:
        ped.writelines(lines)
        ped.flush()

def get_sample_id_from_bam(bam_path):
    """ extracts a sample ID from a BAM file
    
    Args:
        bam_path: path to BAM file
    
    Returns:
        sample ID as string
    """
    
    bam = pysam.AlignmentFile(bam_path)
    read_group = bam.header["RG"]
    
    # select the first read group, and the sample ("SM") value
    sample_id = read_group[0]["SM"]
    
    return sample_id
