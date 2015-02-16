""" script to identify mosaic de novos using denovogear
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os
import subprocess

from src.mosaic_functions import symlink_bam, make_seq_dic_file, \
    make_ped_for_trio, get_sanger_ids, call_mosaic_de_novos
from src.extract_bam import find_bam_on_lustre
from src.mosaic_calling_denovogear import MosaicCalling

PROBANDS_FILE = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/data/probands_without_diagnoses.txt"
# TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bams"
FAMILIES_PED_FILE = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"

def open_families(probands_filename, ped_filename):
    """ find the sample IDs for members of the trios
    """
    
    families = []
    with open(ped_filename) as handle:
        for line in handle:
            if line.startswith("family_id"):
                continue
            
            line = line.strip().split("\t")
            
            # ignore nonproband lines
            paternal_id = line[2]
            if paternal_id == "0":
                continue
            
            family = {}
            family["family_id"] = line[0]
            family["child"] = line[1]
            family["mother"] = line[3]
            family["father"] = paternal_id
            sex = line[4]
            
            families.append([family, sex])
    
    return families

def get_repeats():
    
    repeats = []
    with open("samples_for_rerunning.tsv") as handle:
        for line in handle:
            line = line.strip().split("\t")
            sample_id = line[0]
            repeats.append(sample_id)
    
    return repeats

def main():
    
    all_sanger_ids = get_sanger_ids()
    
    families = open_families(PROBANDS_FILE, FAMILIES_PED_FILE)
    repeats = get_repeats()
    
    # sort the families by DDD sample ID, so that we can reliably get specific
    # sets of families to run
    families = sorted(families, key=lambda k: k[0]["child"])
    
    families = [family for family in families if family[0]["child"] in repeats]
    
    # families = families[0:500] # done on long, many failed, figure out chroms to repeat
    # families = families[500:1000] # run on normal
    # families = families[1000:1500] # run on normal
    # families = families[1500:2000] # run on normal
    # families = families[2000:2500] # run on normal
    # families = families[2500:3000] # run on normal
    # families = families[3000:3500] # run on normal
    # families = families[3500:4000] # run on normal
    # families = families[4000:] # run on normal
    
    # temp = families[0]
    # family = temp[0]
    # sex = temp[1]
    #
    # child_id = family["child"]
    # mother_id = family["mother"]
    # father_id = family["father"]
    
    # child_bam = find_bam_on_lustre(family["child"], all_sanger_ids)
    # mother_bam = find_bam_on_lustre(family["mother"], all_sanger_ids)
    # father_bam = find_bam_on_lustre(family["father"], all_sanger_ids)
    # outdir = os.path.dirname(child_bam)
    
    # # # call a single region
    # caller = MosaicCalling(child_bam, mother_bam, father_bam, sex, outdir)
    # chrom = "1"
    # start = "1"
    # stop = "50000"
    
    # a = "python mosaic_calling_denovogear.py --proband-bam {0} --mother-bam {1} --father-bam {2} --proband-sex {3} --outdir {4} --chrom {5} --start {6} --stop {7}".format(child_bam, mother_bam, father_bam, sex, outdir, chrom, start, stop)
    # print(a)
    #
    # region = (chrom, start, stop)
    # caller.call_mosaic_de_novos_in_region(region)
    #
    # # submit jobs to the cluster to call de novos
    # call_mosaic_de_novos(family, sex, all_sanger_ids)
    
    for family, sex in families:
        call_mosaic_de_novos(family, sex, all_sanger_ids)
    

if __name__ == '__main__':
    main()
