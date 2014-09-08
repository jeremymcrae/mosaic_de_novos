""" small script to correct certain values for the PL field in VCF files.
"""

import os
import sys
import argparse

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description="Fix PL fields in a VCF file.")
    parser.add_argument("infile", nargs="?", type=argparse.FileType("r"), \
        default=sys.stdin, help="VCF to fix")
    
    args = parser.parse_args()
    
    return args.infile


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
    
    try:
        first_line = line
    except:
        first_line = ""
    
    return header, first_line

def process_vcf_line(line):
    """ Fix the VCF PL fields.
    
    Args:
        line: variant line from a VCF file
    
    Returns:
        nothing, instead writes the line to standard out
    """
    
    line = line.rstrip().split("\t")
    format = line[8]
    format_fields = format.split(":")
    
    pl_pos = format_fields.index("PL")
    
    # pull out the PL field for each sample, and where it has a "." in the PL
    # field, we change it to an extremely low phred-scaled likelihood (255)
    for sample_pos in range(9, len(line)):
        sample = line[sample_pos]
        if sample == ".":
            continue
        
        fields = sample.split(":")
        fields[pl_pos] = fields[pl_pos].replace(".", "255")
        
        line[sample_pos] = ":".join(fields)
    
    line = "\t".join(line) + "\n"
    sys.stdout.write(line)

def fix_merged_vcf(merged):
    """ Fix the VCF PL fields. This step is necessary because when you 
    merged the three VCFs, in some cases there may not be enough information
    to determine genotype likelihoods for genotypes involving an alternative
    allele in all three samples. In such a case Vcftools' merge function 
    substitutes the genotype likelihood in the resulting VCF's PL field 
    with a dot ("."). This is not something DNG can work with, and therefore
    the dots have to be substituted with extremely low phred-scaled genotype 
    likelihoods such as 255. This is what the plFix.pl script does. Please 
    note that this script is quite crude and you should probably understand 
    what it does before using it for something important.
    """
    
    # sometimes we just pass in the standard input, rather than a filename. The
    # standard input comes as an opened file-like interface, but given a 
    # filename, the filename needs to be opened.
    if not isinstance(merged, file):
        f = open(merged, "r")
    else:
        f = merged
    
    header, first_line = get_header(f)
    sys.stdout.writelines(header)
    
    process_vcf_line(first_line)
    for line in f:
        process_vcf_line(line)

def main():
    
    infile = get_options()
    fix_merged_vcf(infile)
    

if __name__ == '__main__':
    main()
