""" basic filtering of mosaic de novos using overlap with the standard 
samtools/denovogear output.
"""

from __future__ import print_function
from __future__ import division

import argparse
import sys

def get_options():
    """
    """
    
    
    parser = argparse.ArgumentParser(description="Fix PL fields in a VCF file.")
    parser.add_argument("--standard", help="denovogear output from standard samtools")
    parser.add_argument("--modified", help="denovogear output from modified samtools")
    
    args = parser.parse_args()
    
    return args.standard, args.modified

def parse_denovogear_variant_line(line):
    """ parses a denovogear variant output line
    
    Args:
        line: line for a variant in a denovogear output fields
    
    Returns:
        a dictionary filled with entries from the denovogear line
    """
    
    variant = {}
    
    fields = line.strip().split(" ")
    
    variant["type"] = fields[0]
    
    read_idx = fields.index("READ_DEPTH")
    map_idx = fields.index("MAPPING_QUALITY")
    
    # grab all the fields prior to the read depth fields, since these are easy
    # unique key, value pairs up until then.
    for i in range(2, read_idx, 2):
        key = fields[i].strip(":")
        variant[key] = fields[i + 1]
    
    # grab the person-specific read depth fields, which are unfortunately named
    # the same as for the mapping quality fields
    variant["READ_DEPTH"] = {}
    for i in range(read_idx + 1, read_idx + 6, 2):
        key = fields[i].strip(":")
        variant["READ_DEPTH"][key] = fields[i + 1]
    
    # grab the person-specific mapping quality fields, which are unfortunately 
    # named the same as for the read depth fields
    variant["MAPPING_QUALITY"] = {}
    for i in range(map_idx + 1, map_idx + 6, 2):
        key = fields[i].strip(":")
        variant["MAPPING_QUALITY"][key] = fields[i + 1]
    
    # might as well hang on to the original line, if we want to simply dump 
    # these for later analysis
    variant["line"] = line
    
    return variant

def read_denovogear_output(path):
    """ opens a denovo gear output file, extracts and parses the variant lines
    
    Args:
        path: path to a denovogear output file
    
    Returns:
        dictionary of parsed variant lines, indexed by (chrom, pos) tuples
    """
    
    variants = {}
    
    with open(path, "r") as f:
        for line in f:
            if line.startswith("DENOVO"):
                var = parse_denovogear_variant_line(line)
                key = (var["ref_name"], var["coor"])
                variants[key] = var
    
    return variants

def get_mosaic_only_de_novos(standard, modified):
    """ finds the variants that only occur in the mosaic analyses
    """
    
    mod = read_denovogear_output(modified)
    std = read_denovogear_output(standard)
    
    mosaic_keys = set(mod) - set(std)
    mosaic_only = {key: mod[key] for key in mosaic_keys}
    
    return mosaic_only
        
def main():
    
    standard, modified = get_options()
    
    mosaic = get_mosaic_only_de_novos(standard, modified)
    
    for key in mosaic:
        var = mosaic[key]
        depth = var["READ_DEPTH"]
        qual = var["MAPPING_QUALITY"]
        out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(var["ref_name"], var["coor"], var["pp_dnm"], depth["child"], depth["mom"], depth["dad"], qual["child"], qual["mom"], qual["dad"])
        sys.stdout.write(out)

if __name__ == '__main__':
    main()

    