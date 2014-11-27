""" a class to construct and run the mosaic de novo calling commands for a trio
"""

import os
import argparse
import subprocess
import datetime
import logging

from mosaic_functions import get_sample_id_from_bam, make_corrected_vcf_header, \
    make_ped_for_trio, symlink_bam, make_seq_dic_file

logging.basicConfig(filename='mosaic_calling.log',level=logging.DEBUG)

def get_options():
    """ parse the command line options for the script
    """
    
    parser = argparse.ArgumentParser(description="Run mosaic denovo calling.")
    parser.add_argument("--proband-bam", help="BAM file for proband")
    parser.add_argument("--mother-bam", help="BAM file for mother")
    parser.add_argument("--father-bam", help="BAM File for father")
    parser.add_argument("--proband-sex", help="Gender of proband")
    
    # and define the region of the genome to call
    parser.add_argument("--chrom", help="Chromosome to find denovos in")
    parser.add_argument("--start", help="Region of chromosome to start examining for de novos")
    parser.add_argument("--stop", help="Region of chromosome to stop examining for de novos")
    
    args = parser.parse_args()
    
    region = (args.chrom, args.start, args.stop)
    
    return args.proband_bam, args.mother_bam, args.father_bam, \
        args.proband_sex, region

class MosaicCalling(object):
    """ class to construct and run the mosaic de novo calling commands for a trio
    """
    
    hgi = "/software/hgi/pkglocal"
    
    denovogear = "/nfs/users/nfs_s/sa9/scripts/denovogear/denovogear-0.5/build/src/denovogear"
    reference = "/software/ddd/resources/v1.2/hs37d5.fasta"
    standard_samtools = os.path.join(hgi, "samtools-0.1.19", "bin", "samtools")
    modified_samtools = "/nfs/team29/aw15/samtoolsMod/samtools"
    old_bcftools = "/software/vertres/bin-external/bcftools-0.1.18"
    new_bcftools = os.path.join(hgi, "bcftools-1.1", "bin", "bcftools") # version 1.0+ is necessary for the reheader command
    pl_fixer = os.path.join(os.path.dirname(__file__), "fix_pl_field.py")
    
    female_codes = ["F", "Female", "female", "2"]
    male_codes = ["M", "Male", "male", "1"]
    
    dic_path = "seq_dic.txt"
    make_seq_dic_file(dic_path)
         
    def __init__(self, child_bam, mother_bam, father_bam, sex):
        """ initiates the class with the BAM paths etc
        
        Args:
            child_bam: path to proband's BAM file.
            mother_bam: path to mother's BAM file.
            father_bam: path to father's BAM file.
            sex: sex of the probanddic_path: path to sequence dictionary file
        """
        
        if sex in self.female_codes:
            self.proband_sex = "2"
        elif sex in self.male_codes:
            self.proband_sex = "1"
        else:
            raise ValueError("unknown gender: " + sex)
        
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
        
        # make header files that provide more info than standard samtools output
        # which is required by a recent bcftools version.
        make_corrected_vcf_header(self.child_bam)
        make_corrected_vcf_header(self.mother_bam)
        make_corrected_vcf_header(self.father_bam)
        make_corrected_vcf_header(self.new_bam)
    
    def call_mosaic_de_novos_in_region(self, region):
        """ call the rest of the functions in this class, in the correct order, 
        and on the correct files
        
        Args:
            region: tuple of (chrom, start nucleotide, stop nucleotide) strings
                that define the region of the genome to examine for mosaic de novos
        
        Returns:
            nothing
        """
        
        # run samtools, with the modified samtools used for the child
        self.samtools(self.child_bam, region, modified=True)
        self.samtools(self.mother_bam, region)
        self.samtools(self.father_bam, region)
        self.samtools(self.new_bam, region)
        
        # prepare a BCF file for denovogear, then run denovogear on that
        self.run_denovogear(self.child_bam, self.mother_bam, self.father_bam, region, "standard")
        self.run_denovogear(self.new_bam, self.mother_bam, self.father_bam, region, "modified")
        
        # and tidy up the VCFs that were produced (but leave the BCFs?)
        self.remove_vcf(self.child_bam, region)
        self.remove_vcf(self.mother_bam, region)
        self.remove_vcf(self.father_bam, region)
        self.remove_vcf(self.new_bam, region)
    
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
        region_path = ".{0}.{1}-{2}".format(*region)
        logging.info("calling genotypes for {0} in {1}".format(region_path, bam))
        
        assert os.path.exists(bam)
        assert os.path.exists(self.reference)
        
        sam_cmd = [self.standard_samtools, "mpileup", "-r", region_id, "-gDf", \
            self.reference, bam]
        if modified:
            sam_cmd[0] = self.modified_samtools
            sam_cmd.append("-p0.25")
        
        samtools = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, \
            stderr=open(os.devnull, "w"))
        
        vcf_convert = self.convert_to_vcf(bam, region, samtools)
    
    def convert_to_vcf(self, bam, region, samtools):
        """ convert samtools mpileup output to VCF
        
        Args:
            bam: path to bam file
            region: tuple of (chrom, start, stop) strings
            samtools: Popen command from running samtools, contains samtools 
                output in a pipe.
        
        Returns:
            nothing
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        logging.info("converting samtools output to vcf for {0} in {1}".format(region_path, bam))
        
        # set the path to the output VCF
        temp_vcf = os.path.splitext(bam)[0] + region_path + ".temp.vcf.gz"
        vcf = os.path.splitext(bam)[0] + region_path + ".vcf.gz"
        header = os.path.join(os.path.dirname(bam), "fixed_header.txt")
        
        # convert the samtools output to VCF
        to_vcf = subprocess.Popen([self.old_bcftools, "view", "-"], \
            stdin=samtools.stdout, stdout=subprocess.PIPE)
        
        # compress and tabix the initial VCF output
        bgzip = subprocess.call(["bgzip"], stdin=to_vcf.stdout, stdout=open(temp_vcf, "w"))
        tabix = subprocess.call(["tabix", "-f", "-p", "vcf", temp_vcf])
        
        # fix the header of the initial VCF output
        reheader = subprocess.call([self.new_bcftools, "reheader", "--header", \
            header, temp_vcf], stdout=open(vcf, "w"))
        new_tabix = subprocess.call(["tabix", "-f", "-p", "vcf", vcf])
        
        # and remove the temporary files
        os.remove(temp_vcf)
        os.remove(temp_vcf + ".tbi")
    
    def run_denovogear(self, child, mother, father, region, modify_string):
        """ merge the VCFs from the trio members into a multi-sample VCF
        
        Args:
            child: path to proband's BAM file.
            mother: path to mother's BAM file.
            father: path to father's BAM file.
            region: tuple of (chrom, start, stop) strings
            modify_string: either "standard" or "modified" to indicate whether 
                the proband's data was generated using the standard or modified
                samtools (the modified samtools allows mosaic calling).
        
        Returns:
            nothing
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        child_dir = os.path.dirname(child)
        basename = os.path.basename(child)
        sample_id = basename.split(".")[0]
        
        bcf = os.path.join(child_dir, sample_id + ".merged" + region_path + "." + modify_string + ".bcf")
        dnm = os.path.join(child_dir, sample_id + ".denovogear" + region_path + "." + modify_string + ".dnm")
        
        self.make_bcf_command(child, mother, father, bcf, region)
        self.construct_dng_command(bcf, dnm, region)
    
    def make_bcf_command(self, child, mother, father, bcf, region):
        """ we need to generate BCF for denovogear to work on
        
        Args:
            child: path to proband's BAM file.
            mother: path to mother's BAM file.
            father: path to father's BAM file.
            bcf: path to output BCF (binary call format) data to.
            region: tuple of (chrom, start, stop) strings
        
        Returns:
            nothing
        """
        
        region_path = "{0}.{1}-{2}".format(*region)
        logging.info("preparing bcf for {0} in {1}".format(region_path, child))
        
        # set the paths to the individual VCFs
        child = "{0}.{1}.vcf.gz".format(os.path.splitext(child)[0], region_path)
        mother = "{0}.{1}.vcf.gz".format(os.path.splitext(mother)[0], region_path)
        father = "{0}.{1}.vcf.gz".format(os.path.splitext(father)[0], region_path)
        
        # merge the trio VCFs
        merge = subprocess.Popen([self.new_bcftools, "merge", father, mother, \
            child], stdout=subprocess.PIPE)
        
        # fix the PL field
        pl_fix = subprocess.Popen(["python", self.pl_fixer], stdin=merge.stdout, \
            stdout=subprocess.PIPE)
        
        # generate a BCF for denovogear
        bcf_convert = subprocess.call([self.old_bcftools, "view", "-D", \
            self.dic_path, "-Sb", "-"], stdin=pl_fix.stdout, stdout=open(bcf, "w"))
    
    def construct_dng_command(self, bcf, dnm, region):
        """ Use denovogear to call de novo mutations in a BCF. 
        
        In this particular case, I use denovogear v0.5. As always with DNG, the 
        appropriate PED file (family.ped) has to be in the working directory.
        
        Args:
            bcf: path to BCF data.
            dnm: path to write denovogear results to.
            region: tuple of (chrom, start, stop) strings
        
        Returns:
            nothing
        """
        
        region_path = "{0}.{1}-{2}".format(*region)
        logging.info("running denovogear for {0} in {1}".format(region_path, child))
        
        dng_chr_type = "auto"
        if region[0] == "X" and self.proband_sex in self.male_codes:
            dng_chr_type = "XS"
        elif region[0] == "X" and self.proband_sex in self.female_codes:
            dng_chr_type = "XD"
        
        dng = subprocess.call([self.denovogear, "dnm", dng_chr_type, "--ped", \
            self.ped_path, "--bcf", bcf], stdout=open(dnm, "w"), \
            stderr=open(os.devnull, "w"))
        
        os.remove(bcf)
    
    def remove_vcf(self, sample_bam, region):
        """ remove the temporary VCF files for the region
        
        This relies on removing a file with the correct VCF filename, generated 
        during the convert_to_vcf() method. This function runs after bcf 
        conversion and denovogear analysis, since we construct two BCFs, which
        means we have to call this after both have completed, rather than 
        cleaning the files up at the end of th function.
        
        Args:
            sample_bam: path to bam file
            region: tuple of (chrom, start, stop) strings
        
        Returns:
            nothing
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the path to the individual VCF
        vcf = os.path.splitext(sample_bam)[0] + region_path + ".vcf.gz"
        
        os.remove(vcf)
        os.remove(vcf + ".tbi")

def main():
    """ runs mosaic calling for a single region of the genome in a single trio
    """
    
    # logging.basicConfig(filename="mosaic_de_novo_calling.log")
    
    proband_bam, mother_bam, father_bam, proband_sex, region = get_options()
    
    caller = MosaicCalling(proband_bam, mother_bam, father_bam, proband_sex)
    
    caller.call_mosaic_de_novos_in_region(region)
    

if __name__ == '__main__':
    main()

