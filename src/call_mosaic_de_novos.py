""" script to identify mosaic de novos using denovogear
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os
import subprocess

from src.mosaic_functions import symlink_bam, make_seq_dic_file, \
    make_ped_for_trio, find_bam_path, call_mosaic_de_novos

PROBANDS_FILE = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/data/probands_without_diagnoses.txt"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bams"
FAMILIES_PED_FILE = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/family_relationships.shared.txt"
BAM_EXTRACTOR = "/nfs/users/nfs_j/jm33/apps/VICAR/python/extract_bam.py"

def open_probands_list(filename):
    """ get a list of sample IDs from a file
    """
    
    probands = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            split_line = line.strip().split("\t")
            proband = split_line[1]
            sex = split_line[2]
            probands[proband] = sex
    
    return probands

def open_families(probands_filename, ped_filename):
    """ find the sample IDs for members of the trios
    """
    
    probands = open_probands_list(probands_filename)
    
    families = []
    with open(ped_filename) as f:
        for line in f:
            line = line.strip().split("\t")
            fam_id = line[0]
            individual_id = line[1]
            paternal_id = line[2]
            maternal_id = line[3]
            
            if individual_id in probands:
                family = {}
                family["child"] = individual_id
                family["mother"] = maternal_id
                family["father"] = paternal_id
                sex = probands[individual_id]
                
                families.append((family, sex))
    
    return families

def extract_bams(families, bams_dir):
    """ make sure we have bam files available for all the members of the families
    """
    
    log_dir = "bam_extraction_logs"
    
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    
    for family, sex in families:
        for individual in family:
            sample_id = family[individual]
            bsub_preamble = ["bsub", "-q", "normal", "-o", os.path.join(log_dir, sample_id + ".bjob_output.txt")]
            extract_command = ["python", BAM_EXTRACTOR, "--sample-id", sample_id, "--dir", bams_dir]
            
            subprocess.call(bsub_preamble + extract_command)

def main():
    
    families = open_families(PROBANDS_FILE, FAMILIES_PED_FILE)
    # extract_bams(families, TEMP_DIR)
    
    dic_path = "seq_dic.txt"
    make_seq_dic_file(dic_path)
    
    
    temp = families[0]
    family = temp[0]
    sex = temp[1]
    
    child_id = family["child"]
    mother_id = family["mother"]
    father_id = family["father"]
    
    child_bam = find_bam_path(child_id, TEMP_DIR)
    mother_bam = find_bam_path(mother_id, TEMP_DIR)
    father_bam = find_bam_path(father_id, TEMP_DIR)
    
    new_bam = child_bam[:-3] + "standard_samtools.bam"
    new_child_bam = symlink_bam(child_bam, new_bam)
    
    # make sure there is a ped file available for the trio
    ped_path = os.path.join(os.path.dirname(child_bam), child_id + ".ped")
    make_ped_for_trio(child_bam, mother_bam, father_bam, sex, ped_path)
    
    # caller = MosaicCalling(child_id, mother_id, father_id, child_bam, mother_bam, father_bam, new_bam, sex, dic_path, ped_path)
    
    # chrom = "1"
    # start = "1"
    # # start = "4990000"
    # stop = "5000000"
    # region = (chrom, start, stop)
    # caller.call_mosaic_de_novos_in_region(region)
    
    
    # submit jobs to the cluster to call de novos
    call_mosaic_de_novos(child_bam, mother_bam, father_bam, new_bam, sex, dic_path, ped_path)
    
    # for family, sex in families:
        
    #     child_id = family["child"]
    #     mother_id = family["mother"]
    #     father_id = family["father"]
        
    #     child_bam = find_bam_path(child_id, TEMP_DIR)
    #     mother_bam = find_bam_path(mother_id, TEMP_DIR)
    #     father_bam = find_bam_path(father_id, TEMP_DIR)
        
    #     new_bam = child_bam[:-3] + "standard_samtools.bam"
    #     new_child_bam = symlink_bam(child_bam, new_bam)
        
    #     # make sure there is a ped file available for the trio
    #     ped_path = os.path.join(os.path.dirname(child_bam), child_id + ".ped")
    #     make_ped_for_trio(child_bam, mother_bam, father_bam, sex, ped_path)
        
    #     caller = MosaicCalling(child_id, mother_id, father_id, child_bam, mother_bam, father_bam, new_bam, sex, dic_path, ped_path)
    #     caller.call_mosaic_de_novos()
    

if __name__ == '__main__':
    main()

