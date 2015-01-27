""" basic filtering of mosaic de novos using overlap with the standard 
samtools/denovogear output.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import sys

import pysam
from scipy.stats import poisson
import tabulate

from src.filtering.parse_denovogear import ParseDenovogear

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description="Filter denovogear calls " + \
        "from calling mosaic SNVs.")
    parser.add_argument("--standard", help="denovogear output from standard samtools")
    parser.add_argument("--modified", help="denovogear output from modified samtools")
    parser.add_argument("--proband-bam", required=True, help="bam file for the proband.")
    parser.add_argument("--mother-bam", required=True, help="bam file for the mother.")
    parser.add_argument("--father-bam", required=True, help="bam file for the father.")
    
    args = parser.parse_args()
    
    return args.standard, args.modified, args.proband_bam, args.mother_bam, args.father_bam

def get_mosaic_only_de_novos(standard, modified):
    """ finds the variants that only occur in the mosaic analyses
    """
    
    modified = ParseDenovogear(modified)
    standard = ParseDenovogear(standard)
    
    mosaic_only = modified.get_subset(modified - standard)
    
    return mosaic_only

def count_bases(bam, chrom, pos, max_coverage=1e10, min_qual=0):
    """ counts reads with different base calls at a chrom position
    
    Args:
        bam: pysam bam object
        chrom: chromosome to use eg "chr1" or "1" depending on how the BAM is 
            set up (specifically, an ID found in the BAMs sequence dictionary).
        pos: base position to count bases at.
        max_coverage: maximum coverage at which we stop tallying the bases
        min_qual: minimum quality score required to include the base call
    
    Returns:
        dictionary of read counts indexed by base calls
    """
    
    assert type(pos) == int
    
    bases = {"A": 0, "G": 0, "C": 0, "T": 0, "N": 0}
    
    # count each base at the required site
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos):
        if pileupcolumn.pos != pos - 1:
            continue
        
        for read in pileupcolumn.pileups:
            if read.alignment.is_duplicate: # ignore duplicate reads
                print("duplicate")
                continue
            elif read.alignment.is_qcfail:
                print("FAIL")
                continue
            # Only use alignments with cigar strings of only matches (no soft
            # clipped bases or indels)
            elif read.alignment.cigar[0][0] != 0: 
                continue
            # don't check reads after a certain coverage (typically ~300X)
            elif sum(bases.values()) > max_coverage * 5:
                break
            
            # convert the quality score to integer
            qual = read.alignment.query_qualities[read.query_position]
            if IS_PYTHON2: 
                qual = ord(qual)
            
            if qual < min_qual: # ignore low qual reads
                continue
            
            # get the base call as a string
            base = read.alignment.query_sequence[read.query_position]
            bases[base] += 1
    
    return bases

def examine_variants(mosaic, child_bam, mom_bam, dad_bam):
    
    table = []
    
    variants = mosaic.get_variants()
    
    for (chrom, pos) in sorted(variants):
        # this follows the probability model set out in: 
        # Illumina Inc. Illumina Technical Note: Somatic Variant Caller. (2014). 
        # at <http://res.illumina.com/documents/products/technotes/technote_somatic_variant_caller.pdf>
        
        min_qual = 20
        bases = count_bases(child_bam, chrom, int(pos), max_coverage=1000, min_qual=min_qual)
        mom_bases = count_bases(mom_bam, chrom, int(pos), max_coverage=1000, min_qual=min_qual)
        dad_bases = count_bases(dad_bam, chrom, int(pos), max_coverage=1000, min_qual=min_qual)
        
        # the lambda for the poisson distribution is the proportion of off target
        # bases that we expect, given the read depth and the min_qual
        mu = (10 ** (-(abs(min_qual))/10)) * sum(bases.values())
        
        key = (chrom, pos)
        ref_base = variants[key]["ref_base"]
        alt_bases = set(variants[key]["ALT"].split(","))
        
        try:
            alt_bases.remove("X")
        except KeyError:
            continue
        
        pp_dnm = float(variants[key]["pp_dnm"])
        dng_depth = variants[key]["READ_DEPTH"]["child"]
        child_depth = sum(bases.values())
        mom_depth = sum(mom_bases.values())
        dad_depth = sum(dad_bases.values())
        
        print(child_depth, mom_depth, dad_depth)
        
        # currently hard code some depth filters (maybe swap this to based on the 
        # global coverage)
        if child_depth > 200 or child_depth < 20:
            continue
        
        for alt_base in alt_bases:
            alt_reads = bases[alt_base]
            mom_prp = mom_bases[alt_base]/mom_depth
            dad_prp = dad_bases[alt_base]/dad_depth
            
            de_novo_reads = alt_reads
            if alt_reads > child_depth - alt_reads:
                de_novo_reads = child_depth - alt_reads
                mom_prp = 1 - mom_prp
                dad_prp = 1 - dad_prp
            
            poisson_p = 1 - poisson.cdf(de_novo_reads - 1, mu)
            
            if poisson_p > 0.0001 or pp_dnm < 0.9:
                continue
            # exclude variants with allel frequencies too close to 0.5, since
            # these should have been screened by the standard de novo variant
            # calling pipeline
            if abs((de_novo_reads/child_depth )- 0.5) < 0.05:
                continue
            
            line = [chrom, pos, ref_base, alt_base, poisson_p, pp_dnm, mu, de_novo_reads, mom_prp, dad_prp, dng_depth, child_depth]
            table.append(line)

    header = ["chrom", "pos", "ref", "alt", "poisson_p", "pp_dnm", "mu", \
        "de_novo_reads", "mom_prp", "dad_prp", "dng_depth", "child_depth"]
    
    print(tabulate.tabulate(table, headers = header))

def main():
    
    standard, modified, child_bam_path, mom_bam_path, dad_bam_path = get_options()
    
    child_bam = pysam.AlignmentFile(child_bam_path)
    mom_bam = pysam.AlignmentFile(mom_bam_path)
    dad_bam = pysam.AlignmentFile(dad_bam_path)

    # mosaic = get_mosaic_only_de_novos(standard, modified)
    mosaic = ParseDenovogear(modified)
    
    examine_variants(mosaic, child_bam, mom_bam, dad_bam)
    
    sys.stdout.write("chrom\tpos\tpp_dnm\tchild_read_depth\tmom_read_depth\tdad_read_depth\tchild_qual\tmom_qual\tdad_qual\n")
    
    variants = mosaic.get_variants()
    for key in variants:
        var = variants[key]
        depth = var["READ_DEPTH"]
        qual = var["MAPPING_QUALITY"]
        out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(var["ref_name"], var["coor"], var["pp_dnm"], depth["child"], depth["mom"], depth["dad"], qual["child"], qual["mom"], qual["dad"])
        sys.stdout.write(out)

if __name__ == '__main__':
    main()

    