""" a class to construct and run the mosaic de novo calling commands for a trio
"""

import os
import argparse
import subprocess

from mosaic_functions import get_sample_id_from_bam, make_corrected_vcf_header


def get_options():
    """ parse the command line options for the script
    """
    
    parser = argparse.ArgumentParser(description="Run mosaic denovo calling.")
    parser.add_argument("--proband-bam", help="BAM file for proband")
    parser.add_argument("--mother-bam", help="BAM file for mother")
    parser.add_argument("--father-bam", help="BAM File for father")
    parser.add_argument("--alt-child-bam", help="Alternate BAM file for proband (symlinked from proband BAM)")
    parser.add_argument("--proband-sex", help="Gender of proband")
    parser.add_argument("--sequence-dict", help="path to sequence dictionary")
    parser.add_argument("--ped", help="Path to pedigree file showing trio relationships")
    
    parser.add_argument("--chrom", help="Chromosome to find denovos in")
    parser.add_argument("--start", help="Region of chromosome to start examining for de novos")
    parser.add_argument("--stop", help="Region of chromosome to stop examining for de novos")
    
    args = parser.parse_args()
    
    region = (args.chrom, args.start, args.stop)
    
    return args.proband_bam, args.mother_bam, args.father_bam, \
        args.alt_child_bam, args.proband_sex, args.sequence_dict, args.ped, region

class MosaicCalling(object):
    """ class to construct and run the mosaic de novo calling commands for a trio
    """
    
    denovogear = "/nfs/users/nfs_s/sa9/scripts/denovogear/denovogear-0.5/build/src/denovogear"
    reference = "/software/ddd/resources/v1.2/hs37d5.fasta"
    standard_samtools = "/software/vertres/bin-external/samtools-0.1.18"
    modified_samtools = "/nfs/team29/aw15/samtoolsMod/samtools"
    old_bcftools = "/software/vertres/bin-external/bcftools-0.1.18"
    new_bcftools = "/nfs/users/nfs_j/jm33/apps/bcftools/bcftools" # version 1.0 is necessary for the reheader command
    pl_fixer = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/src/fix_pl_field.py"
    
    female_codes = ["F", "Female", "female", "2"]
    male_codes = ["M", "Male", "male", "1"]
         
    def __init__(self, child_bam, mother_bam, father_bam, new_child_bam, sex, \
        dic_path, ped_path):
        
        if sex in self.female_codes:
            self.proband_sex = "2"
        elif sex in self.male_codes:
            self.proband_sex = "1"
        else:
            raise ValueError("unknown gender: " + sex)
        
        self.dic_path = dic_path
        self.ped_path = ped_path
        
        # find the sample BAMs
        self.child_bam = child_bam
        self.mother_bam = mother_bam
        self.father_bam = father_bam
        
        # make a new bam for the child for using with the standard samtools, so
        # it has a different filename
        self.new_bam = new_child_bam
        
        # make header files that provide more info than standard samtools output
        # which is required by a recent bcftools version.
        make_corrected_vcf_header(self.child_bam)
        make_corrected_vcf_header(self.mother_bam)
        make_corrected_vcf_header(self.father_bam)
        make_corrected_vcf_header(self.new_bam)
    
    def call_mosaic_de_novos_in_region(self, region):
        """ call the rest of the functions in this class, in the correct order, 
        and on the correct files
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
            bcf: path to bcf file (usually located in same folder as bam)
            region: tuple of chrom, start and stop
        
        Returns:
            bsub job ID for the genotype calling cluster job
        """
        
        print(bam)
        
        region_id = "{0}:{1}-{2}".format(*region)
        region_path = ".{0}.{1}-{2}".format(*region)
        
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
            samtools: Popen command from running samtools, contains the output
                as a pipe.
        
        Returns:
            nothing
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the path to the output VCF
        temp_vcf = os.path.splitext(bam)[0] + region_path + ".temp.vcf.gz"
        vcf = os.path.splitext(bam)[0] + region_path + ".vcf.gz"
        header = os.path.dirname(bam) + "fixed_header.txt"
        
        # convert the samtools output to VCF
        to_vcf = subprocess.Popen([self.old_bcftools, "view", "-"], \
            stdin=samtools.stdout, \
            stdout=subprocess.PIPE)
        
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
        """
        
        print(modify_string)
        
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
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the paths to the individual VCFs
        child = os.path.splitext(child)[0] + region_path + ".vcf.gz"
        mother = os.path.splitext(mother)[0] + region_path + ".vcf.gz"
        father = os.path.splitext(father)[0] + region_path + ".vcf.gz"
        
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
        """ Use DNG to call de novo mutations. In this particular case, I use 
        v0.5. As always with DNG, the appropriate PED file (family.ped) has to 
        be in the working directory.
        """
        
        dng_chr_type = "auto"
        if region[0] == "X" and self.proband_sex in self.male_codes:
            dng_chr_type = "XS"
        elif region[0] == "X" and self.proband_sex in self.female_codes:
            dng_chr_type = "XD"
        
        dng = subprocess.call([self.denovogear, "dnm", dng_chr_type, "--ped", \
            self.ped_path, "--bcf", bcf], stdout=open(dnm, "w"), \
            stderr=open(os.devnull, "w"))
    
    def remove_vcf(self, sample_bam, region):
        """ remove the temporary VCF files for the region
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the paths to the individual VCF
        vcf = os.path.splitext(sample_bam)[0] + region_path + ".vcf.gz"
        
        os.remove(vcf)
        os.remove(vcf + ".tbi")

def main():
    proband_bam, mother_bam, father_bam, alt_child_bam, proband_sex, \
        sequence_dict, ped, region = get_options()
        
    caller = MosaicCalling(proband_bam, mother_bam, father_bam, alt_child_bam, \
        proband_sex, sequence_dict, ped)
    
    caller.call_mosaic_de_novos_in_region(region)
    

if __name__ == '__main__':
    main()

