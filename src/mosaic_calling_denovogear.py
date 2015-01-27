""" a class to construct and run the mosaic de novo calling commands for a trio
"""

import os
import sys
import argparse
import subprocess
import time
import logging
import tempfile

from mosaic_functions import get_sample_id_from_bam, make_corrected_vcf_header, \
    make_ped_for_trio, symlink_bam, make_seq_dic_file, chrom_lengths

logging.basicConfig(filename='mosaic_calling.log',level=logging.DEBUG)

def get_options():
    """ parse the command line options for the script
    """
    
    parser = argparse.ArgumentParser(description="Run mosaic denovo calling.")
    parser.add_argument("--proband-bam", required=True, help="BAM file for proband")
    parser.add_argument("--mother-bam", required=True, help="BAM file for mother")
    parser.add_argument("--father-bam", required=True, help="BAM File for father")
    parser.add_argument("--proband-sex", required=True, \
        choices=["1", "M", "male", "2", "F", "female"], help="Sex of proband")
    
    # and define the region of the genome to call
    parser.add_argument("--chrom", required=True, help="Chromosome to find denovos in")
    parser.add_argument("--start", help="Region of chromosome to start examining for de novos, omit to process full chromosome")
    parser.add_argument("--stop", help="Region of chromosome to stop examining for de novos, omit to process full chromosome")
    
    parser.add_argument("--outdir", help="Folder to place denovogear results into")
    
    args = parser.parse_args()
    
    region = (args.chrom, args.start, args.stop)
    
    return args.proband_bam, args.mother_bam, args.father_bam, \
        args.proband_sex, region, args.outdir

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
    pl_fixer = os.path.join(os.path.dirname(__file__), "fix_pl_field.py")
    
    # define the genome reference file (perhaps the DDD file isn't available to
    # everyone, in which case swap to reference on my lustre folder).
    reference = "/software/ddd/resources/v1.2/hs37d5.fasta"
    if not os.path.exists(reference):
        reference = "/lustre/scratch113/teams/hurles/users/jm33/hs37d5.fasta"
        
    female_codes = ["f", "female", "2"]
    male_codes = ["m", "male", "1"]
    
    dic_path = "seq_dic.txt"
    make_seq_dic_file(dic_path)
         
    def __init__(self, child_bam, mother_bam, father_bam, sex, output_dir=None):
        """ initiates the class with the BAM paths etc
        
        Args:
            child_bam: path to proband's BAM file.
            mother_bam: path to mother's BAM file.
            father_bam: path to father's BAM file.
            sex: sex of the proband
            output_dir: folder to place the output in, or None
        """
        
        self.proband_sex = sex
        assert self.proband_sex in self.male_codes + self.female_codes
        
        # make sure there is a ped file available for the trio
        self.ped_path = os.path.join(os.path.dirname(child_bam), "family.ped")
        make_ped_for_trio(child_bam, mother_bam, father_bam, sex, self.ped_path)
        
        # find the sample BAMs
        self.child_bam = child_bam
        self.mother_bam = mother_bam
        self.father_bam = father_bam
        
        # make a new bam for the child for using with the standard samtools, so
        # it has a different filename
        self.new_bam = self.child_bam[:-3] + "standard_samtools.bam"
        new_child_bam = symlink_bam(self.child_bam, self.new_bam)
        
        self.output_dir = output_dir
        if self.output_dir is None:
            self.output_dir = os.path.dirname(self.child_bam)
    
    def call_mosaic_de_novos_in_region(self, region):
        """ call the rest of the functions in this class, in the correct order,
        and on the correct files
        
        Args:
            region: tuple of (chrom, start nucleotide, stop nucleotide) strings
                that define the region of the genome to examine for mosaic de novos
        
        Returns:
            nothing
        """
        
        # if we haven't specified a start and end for a chromosomal region, then
        # assume we want to process the entire chromosome
        if region[1] is None:
            region = (region[0], 1, region[2])
        if region[2] is None:
            region = (region[0], region[1], chrom_lengths[region[0]])
        
        # run samtools, with the modified samtools used for the child
        child_vcf = self.samtools(self.child_bam, region, modified=True)
        mother_vcf = self.samtools(self.mother_bam, region)
        father_vcf = self.samtools(self.father_bam, region)
        new_child_vcf = self.samtools(self.new_bam, region)
        
        # prepare a BCF file for denovogear, then run denovogear on that
        self.run_denovogear(child_vcf, mother_vcf, father_vcf, region, "standard")
        self.run_denovogear(new_child_vcf, mother_vcf, father_vcf, region, "modified")
        
        # and tidy up the VCFs that were produced
        self.remove_vcfs([child_vcf, mother_vcf, father_vcf, new_child_vcf])
    
    def samtools(self, bam, region, modified=False):
        """ call genotypes from a BAM file using samtools
        
        Args:
            bam: path to bam filename
            region: tuple of chrom, start and stop
            modified: True/False for whether to use the modified samtools
        
        Returns:
            nothing
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
            sam_cmd.append("-p0.25")
        
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
        vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf.gz")
        temp_vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf.gz")
        header = tempfile.NamedTemporaryFile(mode="w")
        make_corrected_vcf_header(bam, header)
        
        # convert the samtools output to VCF
        to_vcf = subprocess.Popen([self.old_bcftools, "view", "-"], \
            stdin=samtools.stdout, stdout=subprocess.PIPE)
        
        # compress and tabix the initial VCF output
        subprocess.call(["bgzip"], stdin=to_vcf.stdout, stdout=temp_vcf)
        
        # fix the header of the initial VCF output
        subprocess.call([self.new_bcftools, "reheader", "--header", header.name, \
            temp_vcf.name], stdout=vcf)
        subprocess.call(["tabix", "-f", "-p", "vcf", vcf.name])
        
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
            dnm: path to output denovogear data to.
            region: tuple of (chrom, start, stop) strings
        
        Returns:
            nothing
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
            self.dic_path, "-Sb", "/dev/stdin"], stdin=pl_fix.stdout,
            stdout=subprocess.PIPE)
        
        # and run de novogear on the output
        dnm = os.path.join(self.output_dir, "{0}.denovogear.{1}.{2}.dnm".format(child_id, region_path, modify))
        subprocess.call([self.denovogear, "dnm", dng_chr_type, "--ped", \
            self.ped_path, "--bcf", "/dev/stdin"], stdin=bcf.stdout,  \
            stdout=open(dnm, "w"), stderr=open(os.devnull, "w"))
    
    def remove_vcfs(self, vcfs):
        """ remove the temporary VCF files for the region
        
        This relies on removing files with the correct VCF filename, generated
        during the convert_to_vcf() method. This function runs after bcf
        conversion and denovogear analysis, since we construct two BCFs, which
        means we have to call this after both have completed, rather than
        cleaning the files up at the end of th function.
        
        Args:
            vcfs: list of file handles for vcf files
        """
        
        for vcf in vcfs:
            vcf.close()
            os.remove(vcf.name + ".tbi")

def main():
    """ runs mosaic calling for a single region of the genome in a single trio
    """
    
    proband_bam, mother_bam, father_bam, proband_sex, region, outdir = get_options()
    
    caller = MosaicCalling(proband_bam, mother_bam, father_bam, proband_sex, outdir)
    
    caller.call_mosaic_de_novos_in_region(region)
    

if __name__ == '__main__':
    main()
