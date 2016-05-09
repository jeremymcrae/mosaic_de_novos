""" a bunch of functions to help with mosaic de novo calling
"""

from __future__ import absolute_import

import os
import pysam

from lustre_bams import find_bam_on_lustre
from lsf_jobs import submit_bsub_job, get_random_string

chrom_lengths = {"1": 248956422, "2": 242193529, "3": 198295559,
    "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, \
    "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622,  \
    "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, \
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, \
    "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895}

TEMP_DIR = "/nfs/users/nfs_j/jm33/temp_mosaic"
BAM_EXTRACTOR = "/nfs/users/nfs_j/jm33/apps/VICAR/python/extract_bam.py"
SANGER_IDS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/person_sanger_decipher.txt"
LOG_FILE = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/mosaic_calling.log"

def call_mosaic_de_novos(family, sex, all_sanger_ids=None):
    """ run through all of the chroms
    """
    
    if all_sanger_ids is None:
        all_sanger_ids = get_sanger_ids()
    
    # find the location of the BAMs for the family members
    child_bam = find_bam_on_lustre(family["child"], all_sanger_ids)
    mother_bam = find_bam_on_lustre(family["mother"], all_sanger_ids)
    father_bam = find_bam_on_lustre(family["father"], all_sanger_ids)
    
    # make sure there is a output folder available for the proband
    child_id = get_sample_id_from_bam(child_bam)
    outdir = os.path.join(TEMP_DIR, child_id)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    commands = []
    for chrom in sorted(chrom_lengths):
        job_id = "mosaic_calling_{0}-chr{1}".format(child_id, chrom)
        
        command = ["python3", "src/mosaic_calling_denovogear.py", \
            "--proband-bam", child_bam, \
            "--mother-bam", mother_bam, \
            "--father-bam", father_bam, \
            "--proband-sex", sex, \
            "--chrom", chrom, \
            "--outdir", outdir]
        
        command = " ".join(command)
        commands.append(command)
    
    commands += get_merge_merge_denovogear(outdir, child_id)
    
    command = ["\n".join(commands)]
    submit_bsub_job(command, job_id, memory=500, requeue_code=99, queue="long", cpus=2)

def get_merge_merge_denovogear(folder, child_id):
    """ merge all the denovogear output files, following denovogear calling
    """
    
    # merge all the denovogear output for the modified samtools
    merge = ["python3", "src/filtering/merge_denovogear.py", \
        "--folder", folder, \
        "--remove-files", \
        "--pattern", "{0}*modified".format(child_id), \
        ">", os.path.join(folder, "{0}.denovogear.modified.dnm".format(child_id))]
    
    command = [" ".join(merge)]
    
    return command

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

def get_sanger_ids():
    """ make a dictionary of sanger IDs matched to DDD stable IDs.
    
    Returns:
        dictionary of sanger IDs lists, keyed by their stable ID.
    """
    
    sanger_ids = {}
    with open(SANGER_IDS_PATH) as handle:
        for line in handle:
            line = line.strip().split("\t")
            stable_id = line[0]
            sanger_id = line[2]
            
            if stable_id not in sanger_ids:
                sanger_ids[stable_id] = []
            
            sanger_ids[stable_id].append(sanger_id)
    
    return sanger_ids

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