""" adjusts read group sample information in a BAM merged from BAMs belonging to
two different individuals. 

This script has an odd behaviour, where it requires the merged BAM file to be
passed in as standard input, due to a quirk of how pysam handles standard input.
Pysam can't use standard input as a python file-like interface, instead it
somehow grabs the standard input internally.

This requires standard input in order to use this in a long series of piped
commands. The merged BAM output is written to standard out.

Usage:

cat INTERMEDIATE_BAM | python fix_read_group.py \
    --subsampled SUBSAMPLED_BAM \
    --alternate ALTERNATE_BAM 
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import argparse
import os
import sys

import pysam

from mosaic_functions import get_sample_id_from_bam


def get_options():
    """ parse the command line options
    """
    
    parser = argparse.ArgumentParser(description="fix the read group sample ID" + \
        " in a merged BAM file.")
    parser.add_argument("--subsampled", help="Path to BAM file that was subsampled from")
    parser.add_argument("--alternate", help="Path to BAM file that the subsampled BAM was merged with")
    
    args = parser.parse_args()
    
    return args.subsampled, args.alternate

def get_bam_header(bam_path, sample_id):
    """ extract a bam header, and amends the sample ID if necessary
    
    Args:
        bam_path: path to a BAM file
        sample_id: sample ID to be added into the sample
    """
    
    bam = pysam.AlignmentFile(bam_path)
    
    for pos in range(len(bam.header["RG"])):
        bam.header["RG"][pos]["SM"] = sample_id
    
    return bam.header.copy()

def get_merged_bam_header(bam1_path, bam2_path, sample_id):
    """ makes a BAM header that is the union of two BAM file headers
    
    Args:
        bam1_path: path to first BAM file
        bam2_path: path to second BAM file
        sample_id: we standardise the sample ID field in the "RG" header to this
    
    Returns:
        BAM header merged from the two BAM files, with the sample ID
        standardised within the header.
    """
    
    # get the headers from the original BAM files, but standardise the sample ID
    header = get_bam_header(bam1_path, sample_id)
    header_2 = get_bam_header(bam2_path, sample_id)
    
    # merge the two headers
    header.update(header_2)
    
    return header

def fix_merged_bam(subsampled_path, alternate_path):
    """ creates a merged BAM file from standard input
    
    This function takes a BAM passed in from standard input, and creates a new
    BAM, with the header modified to incorporate the headers from the source 
    BAMs, except that the sample ID is standardised between the two two source
    BAMs.
    
    Args:
        subsampled_path: path to BAM that was subsampled from
        alternate_path: path to the BAM that the subsampled BAM was merged with
    
    Returns:
        nothing
    """
    
    # figure out a fixed header
    sample_id = get_sample_id_from_bam(subsampled_path)
    merged_header = get_merged_bam_header(subsampled_path, alternate_path, sample_id)
    
    # in addition to the parsed commandline options, we also pass in the 
    # complete merged BAM as standard input, but pysam has difficultly with
    # standard input, so we can't specify it as sys.stdin, instead it needs to 
    # use the pysam syntax of reading a filename of "-".
    input_bam = pysam.AlignmentFile("-", "rb")
    
    # open an output BAM, using the corrected header, and write the standard 
    # input merged BAM to that.
    output_bam = pysam.AlignmentFile("-", "wb", header = merged_header)
    
    for item in input_bam: 
        output_bam.write(item)
    
def main():
    
    subsampled, alternate = get_options()
    
    fix_merged_bam(subsampled, alternate)

if __name__ == '__main__':
    main()
