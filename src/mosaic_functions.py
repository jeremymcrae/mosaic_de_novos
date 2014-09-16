""" a bunch of functions to help with mosaic de novo calling
"""

import os
import random
import pysam
import sys
import subprocess

chrom_lengths = {"1": 248956422, "2": 242193529, "3": 198295559, 
    "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, \
    "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622,  \
    "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, \
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, \
    "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895}

def call_mosaic_de_novos(child_bam, mother_bam, father_bam, sex):
    """ run through all of the chroms, region by region
    """
    
    increment = 50000000
    
    i = 1
    for chrom in chrom_lengths:
        max_length = chrom_lengths[chrom]
        start = 1
        end = 1
        
        while end <= max_length:
            start = end
            end = end + increment
            
            job_id = "{0}:{1}-{2}".format(chrom, start, end) + get_random_string()
            
            command = ["python3", "src/mosaic_calling_denovogear.py", \
                "--proband-bam", child_bam, \
                "--mother-bam", mother_bam, \
                "--father-bam", father_bam, \
                "--proband-sex", sex, \
                "--chrom", chrom, \
                "--start", str(start), \
                "--stop", str(end)]
            
            submit_bsub_job(command, job_id, dependent_id=None, memory=500)
            i += 1
            sys.exit() 

def make_corrected_vcf_header(bam_path):
    """ makes a header file that fixes the lack of explanatory lines
    """
    
    # here's a VCF header that includes the potential VCF output types from
    # samtools
    header = ["##fileformat=VCFv4.1\n",
        "##samtoolsVersion=0.1.18 (r982:295)",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n",
        "##INFO=<ID=I16,Number=.,Type=Integer,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">\n",
        "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n",
        "##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\">\n",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"]
    
    bam_dir = os.path.dirname(bam_path)
    new_path = os.path.join(bam_dir, "fixed_header.txt")
    
    # don't remake the corrected header if the file already exists
    if os.path.exists(new_path):
        return
    
    alt_id = get_sample_id_from_bam(bam_path)
    header[-1] = header[-1] + "\t" + alt_id + "\n"
    
    # and write the final line that includes the alt ID
    output = open(new_path, "w")
    output.writelines(header)
    output.close()

def submit_bsub_job(command, job_id, dependent_id=None, memory=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    job = "-J \"{0}\"".format(job_id)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    dependent = ""
    if dependent_id is not None:
        dependent = "-w \"{0}\"".format(dependent_id)
    
    preamble = ["bsub", job, dependent, "-q", "normal", "-o", "bjob_output.txt", mem]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def symlink_bam(current_path, new_path):
    """ make a new bam for the child for using with the standard samtools,
    so it has a different filename
    
    Args:
        current_path: path to exisitng BAM file
        new_path: path for symlinked file
    
    Returns:
        nothing
    """
    
    # don't relink the BAM if the symlink already exists
    if os.path.exists(new_path):
        return
    
    # allow for if the symlink exists, but doesn't point to a valid path
    if os.path.lexists(new_path):
        os.remove(new_path)
        os.remove(new_path + ".bai")
        
    os.symlink(current_path, new_path)
    os.symlink(current_path + ".bai", new_path + ".bai")

def make_seq_dic_file(dic_path):
    """ make sure we have a contig dictionary file available for bcftools
    
    Args:
        dic_path: path to the sequence dictionary file
        
    Returns:
        nothing
    """
    
    # don't remake the sequence dictionary if the file already exists
    if os.path.exists(dic_path):
        return
    
    chroms = list(range(1, 22)) + ["X", "Y"]
    chroms = [str(x) for x in chroms]
    chroms = "\n".join(chroms) + "\n"
    
    seq_dic = open(dic_path, "w")
    seq_dic.write(chroms)
    seq_dic.close()

def find_bam_path(sample_id, bam_dir):
    """ find the path to the extracted BAM file
    
    Args:
        sample_id: DDD sample ID for individual
        bam_dir: folder wher BAM files are located (within per-individual folders)
    
    Returns:
        path to bam file for sample
    """
    
    sample_dir = os.path.join(bam_dir, sample_id.rstrip(".standard_samtools"))
    sample_bam = os.path.join(sample_dir, sample_id + ".bam")
    
    return sample_bam

def make_ped_for_trio(child, mother, father, sex, ped_path):
    """ make a PED file to define the relationship of samples in a BCF file
    
    Args:
        child: path to a bam file for the child.
        mother: path to a bam file for the mother.
        father: path to a bam file for the father.
        ped_path: path to write the PED file to.
    
    Returns:
        nothing
    """
    
    # only create the file if it doesn't already exist
    if os.path.exists(ped_path):
        return
    
    # get the alternative IDs that will exist in the VCF file
    child_id = get_sample_id_from_bam(child)
    mother_id = get_sample_id_from_bam(mother)
    father_id = get_sample_id_from_bam(father)
    
    # make a ped file
    fam_id = "temp"
    
    # make the lines for the PED file, in PED format.
    child_line = "\t".join([fam_id, child_id, father_id, mother_id, sex, "2"]) + "\n"
    father_line = "\t".join([fam_id, father_id, "0", "0"  , "1", "1"]) + "\n"
    mother_line = "\t".join([fam_id, mother_id, "0", "0", "1", "2"]) + "\n"
    lines = [child_line, father_line, mother_line]
    
    output = open(ped_path, "w")
    output.writelines(lines)
    output.close()

def get_sample_id_from_bam(bam_path):
    """ extracts a sample ID from a BAM file
    
    Args:
        bam_path: path to BAM file
    
    Returns:
        sample ID as string
    """
    
    bam = pysam.Samfile(bam_path)
    read_group = bam.header["RG"]
    
    # select the first read group, and the sample ("SM") value
    sample_id = read_group[0]["SM"]
    
    return sample_id

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string():
    """ make a random string, which we can use for bsub job IDs, so that 
    different jobs do not have the same job IDs.
    """
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    # done't allow the random strings to be equivalent to a number, since 
    # the LSF cluster interprets those differently from letter-containing 
    # strings
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    return hash_string
