""" script to generate a ped file for a trio, given a BCF file for the trio.
"""

import os
import argparse
import subprocess

def get_header(f):
    """ get the header lines from a VCF
    
    Args:
        f: handler for a VCF
    
    Returns:
        list of header lines from VCF
    """
    
    header = []
    for line in f:
        if line.startswith("#"):
            header.append(line)
        else:
            break
    
    # and jump back to the start of the file, so other functions can process it
    f.seek(0)
    
    return header

def make_ped_from_trio_bcf(bcf_path, proband_sex, ped_path):
    """ make a PED file to define the samples in a BCF file
    
    Assumes the BCF contains data for members of a trio, in the order of father,
    mother, child. We use the IDs from the last line of the header.
    
    Args:
        bcf_path: path to the BCF file.
        proband_sex: sex of the child (1 = male, 2 = female).
        ped_path: path to write the PED file to.
    """
    
    temp_vcf_path = os.path.splitext(bcf_path)[0] + ".vcf"
    
    # get the first few hundred lines of the BCF, which should be sufficient to
    # capture the full path. Use bcftools to convert the start of the BCF to VCF
    # format.
    vcf_command = ["bcftools", "view", bcf_path, "|", "head", "-n", "300", ">", temp_vcf_path]
    vcf_command = " ".join(vcf_command)
    subprocess.call(vcf_command, shell=True)
    
    # get the sample IDs from the VCF file, which aren't necessarily the
    # same as the sample IDs in the study
    f = open(temp_vcf_path)
    header = get_header(f)
    f.close()
    
    # remove the temporary VCF file
    os.remove(temp_vcf_path)
    
    # get the sample IDs from the header
    header_line = header[-1].rstrip().split("\t")
    sample_ids = header_line[9:]
    father_id, mother_id, child_id = sample_ids
    
    # make a ped file
    fam_id = "temp"
    
    # make the lines for the PED file, in PED format.
    child_line = "\t".join([fam_id, child_id, father_id, mother_id, proband_sex, "2"]) + "\n"
    father_line = "\t".join([fam_id, father_id, "0", "0"  , "1", "1"]) + "\n"
    mother_line = "\t".join([fam_id, mother_id, "0", "0", "1", "2"]) + "\n"
    lines = [child_line, father_line, mother_line]
    
    output = open(ped_path, "w")
    output.writelines(lines)
    output.close()

def get_options():
    """ get the options for the script
    """
    
    parser = argparse.ArgumentParser(description="Make a PED file that represents the individuals in a BCF file.")
    parser.add_argument("--bcf", help="BCF to examine")
    parser.add_argument("--sex", help="sex of the proband")
    parser.add_argument("--ped", help="path to write the ped file to")
    
    args = parser.parse_args()
    
    return args.bcf, args.sex, args.ped


def main():
    bcf_path, proband_sex, ped_path = get_options()
    make_ped_from_trio_bcf(bcf_path, proband_sex, ped_path)
    
if __name__ == '__main__':
    main()

