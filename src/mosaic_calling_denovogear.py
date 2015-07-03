""" a class to construct and run the mosaic de novo calling commands for a trio
"""

from __future__ import division

import os
import sys
import argparse
import subprocess
import time
import logging
import tempfile

from utils import get_sample_id_from_bam, make_corrected_vcf_header, \
    make_ped_for_trio, make_seq_dic_file, chrom_lengths

logging.basicConfig(filename='mosaic_calling.log',level=logging.DEBUG)

def get_options():
    """ parse the command line options for the script
    """
    
    parser = argparse.ArgumentParser(description="Run mosaic denovo calling.")
    parser.add_argument("--proband-bam", required=True, help="BAM file for proband")
    parser.add_argument("--mother-bam", required=True, help="BAM file for mother")
    parser.add_argument("--father-bam", required=True, help="BAM File for father")
    parser.add_argument("--proband-sex", required=True, \
        choices=["1", "M", "m", "Male", "male", "2", "F", "f", "Female", "female"], \
        help="Sex of proband")
    parser.add_argument("--outdir", help="Folder to place denovogear results into")
    
    # and define the region of the genome to call
    parser.add_argument("--chrom", required=True, help="Chromosome to find denovos in")
    parser.add_argument("--start", help="Region of chromosome to start examining for de novos, omit to process full chromosome")
    parser.add_argument("--stop", help="Region of chromosome to stop examining for de novos, omit to process full chromosome")
    
    parser.add_argument("--proportion", type=float, default=0.25, \
        help="Expected proportion of somatic mosaicism")
    
    args = parser.parse_args()
    
    region = (args.chrom, args.start, args.stop)
    
    if args.proportion > 1 or args.proportion < 0:
        sys.exit("error: argument --proportion: the expected proportion of somatic mosaicism must be between 0 and 1.")
    
    return args.proband_bam, args.mother_bam, args.father_bam, \
        args.proband_sex, region, args.outdir, args.proportion

class MosaicCalling(object):
    """ class to construct and run the mosaic de novo calling commands for a trio
    """
    
    hgi = "/software/hgi/pkglocal"
    
    # define all the software tools
    denovogear = "/nfs/users/nfs_s/sa9/scripts/denovogear/denovogear-0.5/build/src/denovogear"
    standard_samtools = "/software/vertres/bin-external/samtools-0.1.18"
    modified_samtools = "/nfs/team29/aw15/samtoolsMod/samtools"
    old_bcftools = "/software/vertres/bin-external/bcftools-0.1.18"
    new_bcftools = os.path.join(hgi, "bcftools-1.1", "bin", "bcftools") # version 1.0+ is necessary for the reheader command
    bgzip = "/software/vertres/bin-external/bgzip"
    tabix = "/software/vertres/bin-external/tabix"
    pl_fixer = os.path.join(os.path.dirname(__file__), "fix_pl_field.py")
    
    # define the genome reference file (perhaps the DDD file isn't available to
    # everyone, in which case swap to reference on my lustre folder).
    reference = "/software/ddd/resources/v1.2/hs37d5.fasta"
    if not os.path.exists(reference):
        reference = "/lustre/scratch113/teams/hurles/users/jm33/hs37d5.fasta"
        
    female_codes = ["f", "female", "2"]
    male_codes = ["m", "male", "1"]
    
    seq_dic = tempfile.NamedTemporaryFile(mode="w")
    make_seq_dic_file(seq_dic)
         
    def __init__(self, child_bam, mother_bam, father_bam, sex, output_dir=None, proportion=0.25):
        """ initiates the class with the BAM paths etc
        
        Args:
            child_bam: path to proband's BAM file.
            mother_bam: path to mother's BAM file.
            father_bam: path to father's BAM file.
            sex: sex of the proband
            output_dir: folder to place the output in, or None
            proportion: expected proportion of somatic mosaic variation.
        """
        
        self.proband_sex = sex.lower()
        assert self.proband_sex in self.male_codes + self.female_codes
        
        # make sure there is a ped file available for the trio
        self.ped = tempfile.NamedTemporaryFile(mode="w")
        make_ped_for_trio(child_bam, mother_bam, father_bam, sex, self.ped)
        
        # find the sample BAMs
        self.child_bam = child_bam
        self.mother_bam = mother_bam
        self.father_bam = father_bam
        
        # catch the temporary vcfs, so that we can delete them later (even if
        # the script crashes or is exited)
        self.temp_vcfs = []
        
        self.output_dir = output_dir
        if self.output_dir is None:
            self.output_dir = os.path.dirname(self.child_bam)
        
        self.proportion = proportion
        assert 0 <= self.proportion <= 1
    
    def call_mosaic_de_novos_in_region(self, region):
        """ call the rest of the functions in this class, in the correct order,
        and on the correct files
        
        Args:
            region: tuple of (chrom, start nucleotide, stop nucleotide) strings
                that define the region of the genome to examine for mosaic de novos
        """
        
        # if we haven't specified a start and end for a chromosomal region, then
        # assume we want to process the entire chromosome
        if region[1] is None:
            region = (region[0], 1, region[2])
        if region[2] is None:
            region = (region[0], region[1], chrom_lengths[region[0]])
        
        try:
            # run samtools, with the modified samtools used for the child
            mother_vcf = self.samtools(self.mother_bam, region)
            father_vcf = self.samtools(self.father_bam, region)
            child_vcf = self.samtools(self.child_bam, region, modified=True)
            
            # prepare a BCF file for denovogear, then run denovogear on that
            self.run_denovogear(child_vcf, mother_vcf, father_vcf, region, "modified")
        finally:
            # ensure we remove the temporary files that were produced
            self.ped.close()
            self.seq_dic.close()
            self.remove_vcfs()
    
    def samtools(self, bam, region, modified=False):
        """ call genotypes from a BAM file using samtools
        
        Args:
            bam: path to bam filename
            region: tuple of chrom, start and stop
            modified: True/False for whether to use the modified samtools
        
        Returns:
            file handle for temporary VCF file
        """
        
        region_id = "{0}:{1}-{2}".format(*region)
        region_path = "{0}.{1}-{2}".format(*region)
        logging.info("{0}\tcalling genotypes for {1} in {2}".format(\
            time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), region_path, bam))
        
        assert os.path.exists(bam)
        assert os.path.exists(self.reference)
        
        bai = bam + ".bai"
        if not os.path.exists(bai):
            subprocess.call([self.standard_samtools, "index", bam])
        
        sam_cmd = [self.standard_samtools, "mpileup", "-r", region_id, "-gDf", \
            self.reference, bam]
        if modified:
            sam_cmd[0] = self.modified_samtools
            sam_cmd.append("-p{0}".format(self.proportion))
        
        samtools = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, \
            stderr=open(os.devnull, "w"))
        
        return self.convert_to_vcf(bam, region, samtools)
    
    def convert_to_vcf(self, bam, region, samtools):
        """ convert samtools mpileup output to VCF
        
        Args:
            bam: path to bam file
            region: tuple of (chrom, start, stop) strings
            samtools: Popen command from running samtools, contains samtools
                mpileup output in a pipe.
        
        Returns:
            handle to VCF file
        """
        
        region_path = "{0}.{1}-{2}".format(*region)
        sample_id = get_sample_id_from_bam(bam)
        logging.info("{0}\tconverting samtools output to vcf for {1} in {2}".format(\
            time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), region_path, bam))
        
        # set the path to the output VCF
        vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf.gz", \
            dir=os.path.dirname(bam), delete=False)
        logging.info("\t{}".format(vcf.name))
        self.temp_vcfs.append(vcf)
        
        temp_vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf.gz")
        header = tempfile.NamedTemporaryFile(mode="w")
        make_corrected_vcf_header(bam, header)
        
        # convert the samtools output to VCF
        to_vcf = subprocess.Popen([self.old_bcftools, "view", "-"], \
            stdin=samtools.stdout, stdout=subprocess.PIPE)
        
        # compress and tabix the initial VCF output
        subprocess.call([self.bgzip], stdin=to_vcf.stdout, stdout=temp_vcf)
        
        # fix the header of the initial VCF output
        subprocess.call([self.new_bcftools, "reheader", "--header", header.name, \
            temp_vcf.name], stdout=vcf)
        subprocess.call([self.tabix, "-f", "-p", "vcf", vcf.name])
        
        # and remove the temp files
        temp_vcf.close()
        header.close()
        
        return vcf
    
    def run_denovogear(self, child, mother, father, region, modify):
        """ we need to merge the VCFs and convert to BCF, then run denovogear
        
        Args:
            child: handle for proband's VCF file.
            mother: handle for mother's VCF file.
            father: handle for father's VCF file.
            region: tuple of (chrom, start, stop) strings.
            modify: either "modified" or "standard" to indicate which samtools
                was used for the proband.
        """
        
        child_id = get_sample_id_from_bam(self.child_bam)
        region_path = "{0}.{1}-{2}".format(*region)
        logging.info("{0}\tmake bcf and run dng for {1} in {2}".format(\
            time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), region_path, child_id))
        
        dng_chr_type = "auto"
        if region[0] == "X" and self.proband_sex in self.male_codes:
            dng_chr_type = "XS"
        elif region[0] == "X" and self.proband_sex in self.female_codes:
            dng_chr_type = "XD"
        
        # merge the trio VCFs
        merge = subprocess.Popen([self.new_bcftools, "merge", father.name, \
             mother.name, child.name], stdout=subprocess.PIPE)
        
        # fix the PL field
        pl_fix = subprocess.Popen(["python", self.pl_fixer], stdin=merge.stdout, \
            stdout=subprocess.PIPE)
        
        # generate a BCF for denovogear
        bcf = subprocess.Popen([self.old_bcftools, "view", "-D", \
            self.seq_dic.name, "-Sb", "/dev/stdin"], stdin=pl_fix.stdout,
            stdout=subprocess.PIPE)
        
        # and run de novogear on the output
        dnm = os.path.join(self.output_dir, "{0}.denovogear.{1}.{2}.dnm".format(child_id, region_path, modify))
        subprocess.call([self.denovogear, "dnm", dng_chr_type, "--ped", \
            self.ped.name, "--bcf", "/dev/stdin"], stdin=bcf.stdout,  \
            stdout=open(dnm, "w"), stderr=open(os.devnull, "w"))
    
    def remove_vcfs(self):
        """ remove the temporary VCF files for the region
        
        This relies on removing files with the correct VCF filename, generated
        during the convert_to_vcf() method. This function runs after bcf
        conversion and denovogear analysis, since we construct two BCFs, which
        means we have to call this after both have completed, rather than
        cleaning the files up at the end of the function.
        
        We can get away with merely closing the vcfs, since they are
        NamedTemporaryFiles, and will be deleted upon closing.
        """
        
        for vcf in self.temp_vcfs:
            vcf.close()
            if os.path.exists(vcf.name + ".tbi"):
                os.remove(vcf.name + ".tbi")

def main():
    """ runs mosaic calling for a single region of the genome in a single trio
    """
    
    proband_bam, mother_bam, father_bam, proband_sex, region, outdir, proportion = get_options()
    
    caller = MosaicCalling(proband_bam, mother_bam, father_bam, proband_sex, outdir, proportion)
    
    try:
        caller.call_mosaic_de_novos_in_region(region)
    except KeyboardInterrupt:
        sys.exit(1)
    

if __name__ == '__main__':
    main()
