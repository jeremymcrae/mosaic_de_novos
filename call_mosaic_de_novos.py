""" script to identify mosaic de novos using denovogear
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os
import subprocess

from src.mosaic_functions import symlink_bam, make_seq_dic_file, \
    make_ped_for_trio, find_bam_path, call_mosaic_de_novos, extract_bams
from src.mosaic_calling_denovogear import MosaicCalling

PROBANDS_FILE = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/data/probands_without_diagnoses.txt"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bams"
FAMILIES_PED_FILE = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/family_relationships.shared.txt"


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
                family["family_id"] = fam_id
                family["child"] = individual_id
                family["mother"] = maternal_id
                family["father"] = paternal_id
                sex = probands[individual_id]
                
                families.append((family, sex))
    
    return families

def main():
    
    families = open_families(PROBANDS_FILE, FAMILIES_PED_FILE)
    # for family, sex in families:
    #     child_id = family["child"]
    #     mother_id = family["mother"]
    #     father_id = family["father"]
    #     for sample_id in [child_id, mother_id, father_id]:
    #          extract_bams(sample_id, TEMP_DIR)
    
    temp = families[0]
    family = temp[0]
    sex = temp[1]
    
    child_id = family["child"]
    mother_id = family["mother"]
    father_id = family["father"]
    
    child_bam = find_bam_path(child_id, TEMP_DIR)
    mother_bam = find_bam_path(mother_id, TEMP_DIR)
    father_bam = find_bam_path(father_id, TEMP_DIR)
    outdir = os.path.dirname(child_bam)
    
    # # call a single region
    caller = MosaicCalling(child_bam, mother_bam, father_bam, sex, outdir)
    chrom = "1"
    start = "1"
    stop = "50000"
    region = (chrom, start, stop)
    caller.call_mosaic_de_novos_in_region(region)
    
    # # submit jobs to the cluster to call de novos
    # call_mosaic_de_novos(family, sex)
    
    # for family, sex in families:
    #     call_mosaic_de_novos(family, sex)
    

if __name__ == '__main__':
    main()

