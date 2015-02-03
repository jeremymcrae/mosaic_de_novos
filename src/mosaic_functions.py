""" a bunch of functions to help with mosaic de novo calling
"""

import os
import random
import pysam
import sys
import subprocess

from extract_bam import get_irods_path_for_participant, extract_bam_from_irods

chrom_lengths = {"1": 248956422, "2": 242193529, "3": 198295559,
    "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, \
    "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622,  \
    "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, \
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, \
    "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895}

TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bams"
BAM_EXTRACTOR = "/nfs/users/nfs_j/jm33/apps/VICAR/python/extract_bam.py"

def call_mosaic_de_novos(family, sex):
    """ run through all of the chroms, region by region
    """
    
    # bams = []
    # bam_ids = []
    # # make sure we have a bam dir
    # for member in ["child", "mother", "father"]:
    #     sample_id = family[member]
    #     bam_paths = get_irods_path_for_participant(sample_id)
    #     if len(bam_paths["lustre"]) == 0:
    #         bam = find_bam_path(child_id, TEMP_DIR)
    #         job_id = extract_bams(sample_id, bam_path=bam)
    #         bam_ids.append(bam_ids)
    #     else:
    #         bam = bam_paths["lustre"][0]
    #     bams.append(bam)
    
    child_bam = find_bam_path(family["child"], TEMP_DIR)
    mother_bam = find_bam_path(family["mother"], TEMP_DIR)
    father_bam = find_bam_path(family["father"], TEMP_DIR)
    
    # check if any of the bam extraction jobs are still running
    
    child_id = get_sample_id_from_bam(child_bam)
    
    job_ids = []
    for chrom in chrom_lengths:
        start = 1
        end = chrom_lengths[chrom]
        
        job_id = "{0}:{1}-{2}".format(chrom, start, end) + get_random_string()
        
        command = ["python3", "src/mosaic_calling_denovogear.py", \
            "--proband-bam", child_bam, \
            "--mother-bam", mother_bam, \
            "--father-bam", father_bam, \
            "--proband-sex", sex, \
            "--chrom", chrom,
            "--start", str(start), \
            "--stop", str(end)]
        
        # submit_bsub_job(command, job_id, dependent_id=bam_ids, memory=500, requeue_code=99)
        submit_bsub_job(command, job_id, memory=500, requeue_code=99)
        job_ids.append(job_id)
    
    merge_denovogear(child_bam, child_id, job_ids)

def merge_denovogear(child_bam, child_id, job_ids):
    """ merge all the denovogear output files, following denovogear calling
    """
    
    # merge all the denovogear output for the standard samtools
    folder = os.path.dirname(os.path.splitext(child_bam)[0])
    job_id = "merge_standard_denovogear" + get_random_string()
    command = ["python3", "src/filtering/merge_denovogear.py", \
        "--folder", folder, \
        "--remove-files", \
        "--pattern", "standard", \
        ">", os.path.join(folder, "{0}.denovogear.standard.dnm".format(child_id))]
    
    submit_bsub_job(command, job_id, dependent_id=job_ids)
    
    # merge all the denovogear output for the modified samtools
    job_id = "merge_modified_denovogear" + get_random_string()
    command = ["python3", "src/filtering/merge_denovogear.py", \
        "--folder", folder, \
        "--remove-files", \
        "--pattern", "modified", \
        ">", os.path.join(folder, "{0}.denovogear.modified.dnm".format(child_id))]
    
    submit_bsub_job(command, job_id, dependent_id=job_ids)

def extract_bams(sample_id, bam_dir=None, bam_path=None):
    """ make sure we have bam files for all the family members
    
    Args:
        family: tuple of (family dict, sex)
        bam_dir: folder to store the BAM file into
        bam_dir: path to store the BAM file at
    """
    
    # only allow one of bam_dir or bam_path to be used
    assert (bam_dir is not None and bam_path is None) or (bam_dir is None and bam_path is not None)
    
    log_dir = "bam_extraction_logs"
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    
    job_id = get_random_string()
    
    log = os.path.join(log_dir, sample_id + ".bjob_output.txt")
    if bam_dir is not None:
        command = ["python", BAM_EXTRACTOR, "--sample-id", sample_id, "--dir", bam_dir]
    elif bam_path is not None:
        command = ["python", BAM_EXTRACTOR, "--sample-id", sample_id, "--path", bam_path]
    
    submit_bsub_job(command, job_id, logfile=log)
    
    return job_id

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

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", "normal", "-o", log, mem]
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

def find_bam_path(sample_id, bam_dir):
    """ find the path to the extracted BAM file
    
    Args:
        sample_id: DDD sample ID for individual
        bam_dir: folder wher BAM files are located (within per-individual folders)
    
    Returns:
        path to bam file for sample
    """
    
    sample_dir = os.path.join(bam_dir, sample_id)
    sample_bam = os.path.join(sample_dir, sample_id + ".bam")
    
    return sample_bam

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
    father_line = "\t".join([fam_id, father_id, "0", "0"  , "1", "1"]) + "\n"
    mother_line = "\t".join([fam_id, mother_id, "0", "0", "1", "2"]) + "\n"
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
