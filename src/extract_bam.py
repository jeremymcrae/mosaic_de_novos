""" functions to extract BAM files for participants from IRODS.

Uses a DDD database tool (adb) to find the IRODS location, and uses an IRODS
tool (iget) to extract the BAM and corresponding BAI file. Requires
authentication via kinit each time this is run, which is made easier by piping
the password from a file (stored as read only for the user) to kinit.

This script can be called as:
    python extract_bam.py --sample-id DDD_SAMPLE_ID --dir DIR_NAME

or used in other scripts with:
from extract_bam import get_irods_path_for_participant, extract_bam_from_irods

"""

from __future__ import print_function
from __future__ import division

import subprocess
import os
import argparse

BAM_PATHS = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/bam_paths.txt"

def get_options():
    """ gets options for the script, if this is being called externally
    """
    
    parser = argparse.ArgumentParser(description="Extract BAM file for a DDD individual.")
    parser.add_argument("--sample-id", required=True, help="individual ID")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dir", help="location to extract bam files to")
    group.add_argument("--path", help="location to extract bam files to")
    
    args = parser.parse_args()
    
    return args.sample_id, args.dir, args.path

def get_irods_path_for_participant(sample_id):
    """ finds the BAM location for a participant from IRODS
    
    Args:
        sample_id: DDD participant ID (eg DDDP100137).
    
    Returns:
        path to the BAM file in IRODs
    """
    
    with open(BAM_PATHS) as handle:
        for line in handle:
            line = line.strip().split("\t")
            person_id = line[0]
            path = line[4]
            if person_id == sample_id:
                return path
    
    return None

def extract_bam_from_irods(irods_bam, output_bam):
    """ extracts BAM and BAI files for a participant from IRODS
    
    Args:
        irods_bam: location of BAM file on IRODS
        output_path: path to extract the BAM file to (lustre recommended).
        attempts: keep track of how many times this function is recursively
            called, so we can quit after too many attempts
    """
    
    # set up the environment to be able to use adb, and access IRODS
    subprocess.call("cat ~/.kinit | kinit", shell=True)
    
    irods_bai = irods_bam + ".bai"
    output_bai = output_bam + ".bai"
    
    # and pull the BAM and BAI files out of IRODs, so that we can work on them
    print("extracting from {} to {}".format(irods_bam, output_bam))
    subprocess.call(["iget", "-v", "-X", "-K", "--retries", "3", irods_bam, output_bam])
    subprocess.call(["iget", "-v", "-X", "-K", "--retries", "3", irods_bai, output_bai])


def main():
    """ extracts the BAM file for a DDD ID, unless it already exists
    
    Args:
        sample_id: DDD participant ID (eg DDDP100137).
        temp_dir: directory in which to store the BAM file.
    
    Returns:
        the path to the BAM file for the participant.
    """
    
    sample_id, output_dir, output_path = get_options()
    
    if output_dir is not None:
        output_dir = os.path.join(output_dir, sample_id)
        output_path = os.path.join(output_dir, sample_id + ".bam")
    else:
        output_dir = os.path.dirname(output_path)
    
    # if the BAM file has been generated already, don't re-request it
    if not os.path.exists(output_path):
        # make sure there's a folder for the BAM (use sample specific folders)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        bam = get_irods_path_for_participant(sample_id)
        extract_bam_from_irods(bam, output_path)
        
    

if __name__ == '__main__':
    main()
