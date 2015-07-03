""" basic filtering of mosaic de novos using overlap with the standard
samtools/denovogear output.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import os
import sys
import subprocess
import json
import tempfile
import itertools

# import pysam
# from scipy.stats import poisson
# import tabulate

from src.filtering.parse_denovogear import ParseDenovogear
from src.utils import get_sanger_ids
# from src.utils import get_sample_id_from_bam

OUTPUT_DIR = "/nfs/users/nfs_j/jm33/mosaic"
DEPTHS_SCRIPT = "/nfs/users/nfs_j/jm33/apps/clinical-filter/scripts/get_trio_depths.py"

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description="Filter denovogear calls " + \
        "from calling mosaic SNVs.")
    parser.add_argument("--denovogear", help="denovogear output")
    
    args = parser.parse_args()
    
    return args.denovogear

def get_ddd_id(denovogear_path):
    """ get the DDD ID for a sample from a denovogear path
    """
    
    basename = os.path.basename(denovogear_path)
    sanger_id = basename.split(".")[0]
    
    sanger_ids = get_sanger_ids()
    
    for key in sanger_ids:
        value = sanger_ids[key]
        
        if sanger_id in value:
            return key

def format_variant(mosaic, key, depths):
    """
    """
    
    var = mosaic[key]
    
    # select the depth entry that is for the required variant
    for depth in depths:
        if depth["chrom"] == var["ref_name"] and str(depth["position"]) == var["coor"]:
            break
    
    # select the correct columns for the variant
    variant_info = [depth["child"]["id"], depth["mom"]["id"], depth["dad"]["id"], var["ID"], \
        depth["chrom"], depth["position"], depth["ref_allele"], depth["alt_allele"], \
        var["type"], var["pp_dnm"]]
    
    # get the forward and reverse depths for the reference and alternate alleles
    # for the child, mother and father
    people = ["child", "mom", "dad"]
    alleles = ["ref", "alt"]
    strands = ["forward", "reverse"]
    combos = itertools.product(people, alleles, strands)
    sample_depths = [ depth[x[0]]["depths"][x[1]][x[2]] for x in combos ]
    
    # and format the lists for writing out
    line = variant_info + sample_depths
    line = [ str(x) for x in line ]
    
    return "\t".join(line) + "\n"

def write_output(output, mosaic, depths):
    
    header = ["child_id", "mom_id", "dad_id", "child_sanger_id", \
        "chrom", "position", "ref", "alt", "type", "PP_DNM", \
        "child_ref_F", "child_ref_R", "child_alt_F", "child_alt_R", \
        "mom_ref_F", "mom_ref_R", "mom_alt_F", "mom_alt_R", \
        "dad_ref_F", "dad_ref_R", "dad_alt_F", "dad_alt_R"]
    
    output.write("\t".join(header) + "\n")
    
    for key in sorted(mosaic):
        line = format_variant(mosaic, key, depths)
        output.write(line)

def main():
    
    denovogear_path = get_options()
    
    proband = get_ddd_id(denovogear_path)
    mosaic = ParseDenovogear(denovogear_path)
    
    # format the variant coordinates, to pass into a script to get depths
    coordinates = []
    for key in mosaic:
        var = mosaic[key]
        entry = {"chrom": var["ref_name"],
            "start": var["coor"],
            "end": var["coor"],
            "ref": var["ref_base"],
            "alt": var["alt_allele"]}
        
        coordinates.append(entry)
    
    # write the variant coordinates to a file
    json_handle = tempfile.NamedTemporaryFile()
    json.dump(coordinates, json_handle, indent=4)
    json_handle.flush()
    
    # get the forward and reverse ref and alt depths for each variant
    depths = subprocess.check_output(["python", DEPTHS_SCRIPT, "--proband", proband, "--json-input", json_handle.name, "--by-strand"])
    depths = json.loads(depths)
    
    output = open(os.path.join(OUTPUT_DIR, "{}_dng_depths.txt".format(proband)), "w")
    write_output(output, mosaic, depths)

if __name__ == '__main__':
    main()
