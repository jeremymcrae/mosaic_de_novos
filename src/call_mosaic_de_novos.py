""" script to identify mosaic de novos using denovogear
"""

from __future__ import print_function
from __future__ import division

import os
import subprocess
import tempfile
import random
import sys
import shutil
import pysam

# from extract_bam import get_irods_path_for_participant, extract_bam_from_irods

PROBANDS_FILE = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/data/probands_without_diagnoses.txt"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bams"
FAMILIES_PED_FILE = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/family_relationships.shared.txt"
BAM_EXTRACTOR = "/nfs/users/nfs_j/jm33/apps/VICAR/python/extract_bam.py"
FIXED_HEADER = os.path.join(TEMP_DIR, "corrected_samtools_vcf_header.txt")
DENOVOGEAR = "/nfs/users/nfs_s/sa9/scripts/denovogear/denovogear-0.5/build/src/denovogear"

def open_probands_list(filename):
    """ get a list of sample IDs from a file
    """
    
    probands = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            split_line = line.strip().split("\t")
            proband = split_line[1]
            sex = split_line[2]
            probands[proband] = sex
    
    return probands

def open_families(probands_filename, ped_filename):
    """ find the sample IDs for members of the trios
    """
    
    probands = open_probands_list(probands_filename)
    
    families = []
    with open(ped_filename) as f:
        for line in f:
            line = line.strip().split("\t")
            fam_id = line[0]
            individual_id = line[1]
            paternal_id = line[2]
            maternal_id = line[3]
            
            if individual_id in probands:
                family = {}
                family["child"] = individual_id
                family["mother"] = maternal_id
                family["father"] = paternal_id
                sex = probands[individual_id]
                
                families.append((family, sex))
    
    return families

def extract_bams(families, bams_dir):
    """ make sure we have bam files available for all the members of the families
    """
    
    log_dir = "bam_extraction_logs"
    
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    
    for family, sex in families:
        for individual in family:
            sample_id = family[individual]
            bsub_preamble = ["bsub", "-q", "normal", "-o", os.path.join(log_dir, sample_id + ".bjob_output.txt")]
            extract_command = ["python", BAM_EXTRACTOR, "--sample-id", sample_id, "--dir", bams_dir]
            
            subprocess.call(bsub_preamble + extract_command)

def get_sample_id_from_bam(bam_path):
    """ extracts a sample ID from a BAM file
    
    Args:
        bam_path: path to BAM file
    
    Returns:
        sample ID as string
    """
    
    bam = pysam.Samfile(bam_path)
    read_group = bam.header["RG"]
    
    # select the first read group, and the sample ("SM") value
    sample_id = read_group[0]["SM"]
    
    return sample_id


class MosaicCalling(object):
    
    reference = "/software/ddd/resources/v1.2/hs37d5.fasta"
    mod_samtools = "/nfs/team29/aw15/samtoolsMod/samtools"
    vcftools_merge = "/software/ddd/external/vcftools/0.1.11/bin/vcf-merge"
    bcftools = "/nfs/users/nfs_j/jm33/apps/bcftools/bcftools"
    pl_fixer = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/src/fix_pl_field.py"
    # ped_maker = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/src/make_ped_from_trio_bcf.py"
    overlap_filter = "/nfs/users/nfs_j/jm33/apps/mosaic_de_novos/src/filter_mosaic_denovogear.py"
    
    female_codes = ["F", "Female", "female", "2"]
    male_codes = ["M", "Male", "male", "1"]
    
    chrom_lengths = {"1": 248956422, "2": 242193529, "3": 198295559, 
        "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, \
        "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622,  \
        "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, \
        "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, \
        "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895, \
        "Y": 57227415}
         
    def __init__(self, family, sex, bams_dir):
        
        self.child_id = family["child"]
        self.mother_id = family["mother"]
        self.father_id = family["father"]
        
        if sex in self.female_codes:
            self.proband_sex = "2"
        elif sex in self.male_codes:
            self.proband_sex = "1"
        else:
            raise ValueError("unknown gender: " + sex)
        
        self.bams_dir = bams_dir
        
        # make sure we have a contig dictionary file available
        self.dic_path = "seq_dic.txt"
        self.make_seq_dic_file()
        
        self.make_corrected_vcf_headers([self.child_id, self.mother_id, self.father_id, self.child_id + ".standard_samtools"])
        
        # find the sample BAMs
        self.child_bam = self.find_bam_path(self.child_id, self.bams_dir)
        self.mother_bam = self.find_bam_path(self.mother_id, self.bams_dir)
        self.father_bam = self.find_bam_path(self.father_id, self.bams_dir)
        
        # make a new bam for the child for using with the standard samtools, so
        # it has a different filename
        self.new_bam = self.child_bam[:-3] + "standard_samtools.bam"
        self.make_new_child_bam()
        
        # make sure there is a ped file available for the trio
        self.ped_path = os.path.join(os.path.dirname(self.child_bam), self.child_id + ".ped")
        self.make_ped_for_trio(self.child_bam, self.mother_bam, self.father_bam, self.ped_path)
    
    def make_seq_dic_file(self):
        """ make sure we have a contig dictionary file available
        """
        
        if not os.path.exists(self.dic_path):
            seq_dic = open(self.dic_path, "w")
            chroms = list(range(1, 22)) + ["X", "Y"]
            chroms = [str(x) for x in chroms]
            chroms = "\n".join(chroms) + "\n"
            seq_dic.write(chroms)
            seq_dic.close()
    
    def make_new_child_bam(self):
        """ make a new bam for the child for using with the standard samtools,
        so it has a different filename
        """
        
        if not os.path.exists(self.new_bam):
            # allow for if the symlink exists, but doesn't point to a valid path
            if os.path.lexists(new_bam):
                os.remove(self.new_bam)
                os.remove(self.new_bam + ".bai")
                
            os.symlink(self.child_bam, self.new_bam)
            os.symlink(self.child_bam + ".bai", self.new_bam + ".bai")
    
    def is_number(self, string):
        """ check whether a string can be converted to a number
        """
        
        try:
            number = float(string)
        except ValueError:
            return False
        
        return True
    
    def get_random_string(self):
        """ make a random string, which we can use for bsub job IDs, so that 
        different jobs do not have the same job IDs.
        """
        
        # set up a random string to associate with the run
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
        
        # done't allow the random strings to be equivalent to a number, since 
        # the LSF cluster interprets those differently from letter-containing 
        # strings
        while self.is_number(hash_string):
            hash_string = "%8x" % random.getrandbits(32)
            hash_string = hash_string.strip()
        
        return hash_string
    
    def find_bam_path(self, sample_id, bam_dir):
        """ find the path to the extracted BAM file
        """
        
        sample_dir = os.path.join(bam_dir, sample_id.rstrip(".standard_samtools"))
        sample_bam = os.path.join(sample_dir, sample_id + ".bam")
        
        return sample_bam
    
    def call_mosaic_de_novos(self):
        """ run through all of the chroms, region by region
        """
        
        increment = 50000000
        
        i = 1
        for chrom in self.chrom_lengths:
            max_length = self.chrom_lengths[chrom]
            start = 1
            end = 1
            
            while end <= max_length:
                start = end
                end = end + increment
                
                region = (chrom, start, end)
                self.call_mosaic_de_novos_in_region(region)
                i += 1
    
    def make_corrected_vcf_headers(self, samples):
        """ makes a header file for each bam that fixes the lack of explanatory lines
        """
        
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        
        for sample_id in samples:
            bam_path = self.find_bam_path(sample_id, self.bams_dir)
            bam_dir = os.path.dirname(bam_path)
            
            alt_id = get_sample_id_from_bam(bam_path)
            header = header_line + "\t" + alt_id + "\n"
            
            new_path = os.path.join(bam_dir, sample_id + ".fixed_header.txt")
            
            shutil.copy(FIXED_HEADER, new_path)
            
            # and write the final line that includes the alt ID
            output = open(new_path, "a")
            output.write(header)
            output.close()
    
    def make_ped_for_trio(self, child, mother, father, ped_path):
        """ make a PED file to define the samples in a BCF file
        
        Assumes the BCF contains data for members of a trio, in the order of father,
        mother, child. We use the IDs from the last line of the header.
        
        Args:
            bcf_path: path to the BCF file.
            ped_path: path to write the PED file to.
        """
        
        # get the alternative IDs that will exist in the VCF file
        child_id = get_sample_id_from_bam(child)
        mother_id = get_sample_id_from_bam(mother)
        father_id = get_sample_id_from_bam(father)
        
        # make a ped file
        fam_id = "temp"
        
        # make the lines for the PED file, in PED format.
        child_line = "\t".join([fam_id, child_id, father_id, mother_id, self.proband_sex, "2"]) + "\n"
        father_line = "\t".join([fam_id, father_id, "0", "0"  , "1", "1"]) + "\n"
        mother_line = "\t".join([fam_id, mother_id, "0", "0", "1", "2"]) + "\n"
        lines = [child_line, father_line, mother_line]
        
        # only create the file if it doesn't already exist
        if not os.path.exists(ped_path):
            output = open(ped_path, "w")
            output.writelines(lines)
            output.close()
    
    def call_mosaic_de_novos_in_region(self, region):
        """ call the rest of the functions in this class, in the correct order, 
        and on the correct files
        """
        
        # run samtools, with the modified samtools used for the child
        child = self.modified_samtools(self.child_bam, region)
        mother = self.standard_samtools(self.mother_bam, region)
        father = self.standard_samtools(self.father_bam, region)
        new_child = self.standard_samtools(self.new_bam, region)
        command = child + [";"] + mother + [";"] + father + [";"] + new_child
        job_id = self.get_random_string() + "_samtools"
        self.submit_bsub_job(command, job_id, memory=300)
        
        # prepare a BCF file for denovogear, then run denovogear on that
        merge_id = self.run_denovogear(self.child_bam, self.new_bam, self.mother_bam, self.father_bam, region, job_id)
    
    def submit_bsub_job(self, command, job_id, dependent_id=None, memory=None):
        """ construct a bsub job submission command
        """
        
        job = "-J \"{0}\"".format(job_id)
        
        mem = ""
        if memory is not None:
            mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
        
        dependent = ""
        if dependent_id is not None:
            dependent = "-w \"{0}\"".format(dependent_id)
        
        preamble = ["bsub", job, dependent, "-q", "normal", "-o", "bjob_output.txt", mem]
        command = ["bash", "-c", "\""] + command + ["\""]
        
        command = " ".join(preamble + command)
        subprocess.call(command, shell=True)
    
    def standard_samtools(self, bam, region):
        """ call genotypes from a BAM file using samtools
        
        Args:
            bam: path to bam filename
            bcf: path to bcf file (usually located in same folder as bam)
            region: tuple of chrom, start and stop
        
        Returns:
            bsub job ID for the genotype calling cluster job
        """
        
        region_id = "{0}:{1}-{2}".format(*region)
        region_path = ".{0}.{1}-{2}".format(*region)
        
        assert os.path.exists(bam)
        assert os.path.exists(self.reference)
        
        bcf = os.path.splitext(bam)[0] + region_path + ".bcf"
        
        job_id = self.get_random_string() + "_samtools"
        samtools = ["samtools", "mpileup", "-r", region_id, "-gDf", self.reference, bam, "|"]
        vcf_convert = self.convert_to_vcf(bam, region)
        
        command = samtools + vcf_convert
        
        return command
    
    def modified_samtools(self, bam, region):
        """ call genotypes from a BAM file using a modified samtools
        
        Args:
            bam: path to bam filename
            bcf: path to bcf file (usually located in same folder as bam)
            region: tuple of chrom, start and stop
        
        Returns:
            bsub job ID for the genotype calling cluster job
        """
        
        assert os.path.exists(self.mod_samtools)
        assert os.path.exists(bam)
        assert os.path.exists(self.reference)
        
        region_id = "{0}:{1}-{2}".format(*region)
        region_path = ".{0}.{1}-{2}".format(*region)
        
        bcf = os.path.splitext(bam)[0] + region_path + ".bcf"
        
        job_id = self.get_random_string() + "_mod_samtools"
        mod_samtools = [self.mod_samtools, "mpileup", "-p0.25", "-r", region_id, "-gDf", self.reference, bam, "|"]
        vcf_convert = self.convert_to_vcf(bam, region)
        
        command = mod_samtools + vcf_convert
        
        return command
    
    def convert_to_vcf(self, bam, region):
        """ convert mpileup-generated BCF files into VCF, then bgzip and tabix
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the path to the output VCF
        vcf = os.path.splitext(bam)[0] + region_path + ".temp.vcf.gz"
        new_vcf = os.path.splitext(bam)[0] + region_path + ".vcf.gz"
        header = os.path.join(os.path.splitext(bam)[0] + ".fixed_header.txt")
        
        # define the unix tools to use
        bcftools = ["bcftools", "view", "-", "|"]
        bgzip = ["bgzip", ">", vcf, ";"]
        tabix = ["tabix", "-f", "-p", "vcf", vcf, ";"]
        reheader = [self.bcftools, "reheader", "--header", header, vcf, ">", new_vcf, ";"]
        new_tabix = ["tabix", "-f", "-p", "vcf", new_vcf, ";"]
        rm_temp = ["rm", vcf, vcf + ".tbi"]
        
        vcf_convert = bcftools +bgzip + tabix + reheader + new_tabix + rm_temp
        
        return vcf_convert
    
    def run_denovogear(self, child, new_child, mother, father, region, samtools_job_id):
        """ merge the VCFs from the trio members into a multi-sample VCF
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        child_dir = os.path.dirname(child)
        merge_mod = os.path.join(child_dir, self.child_id + ".merged" + region_path + ".modified.bcf")
        dnm_mod = os.path.join(child_dir, self.child_id + ".denovogear" + region_path + ".modified.dnm")
        
        merge_std = os.path.join(child_dir, self.child_id + ".merged" + region_path + ".standard.bcf")
        dnm_std = os.path.join(child_dir, self.child_id + ".denovogear" + region_path + ".standard.dnm")
        
        job_id = samtools_job_id.rstrip("samtools") + "denovogear"
        
        bcf_mod = self.make_bcf_command(new_child, mother, father, merge_mod, region)
        bcf_std = self.make_bcf_command(new_child, mother, father, merge_std, region)
        
        dng_mod = self.construct_dng_command(merge_mod, dnm_mod, region)
        dng_std = self.construct_dng_command(merge_std, dnm_std, region)
        
        # set the paths to the individual VCFs
        child_vcf = os.path.splitext(child)[0] + region_path + ".vcf.gz"
        new_child_vcf = os.path.splitext(new_child)[0] + region_path + ".vcf.gz"
        mother_vcf = os.path.splitext(mother)[0] + region_path + ".vcf.gz"
        father_vcf = os.path.splitext(father)[0] + region_path + ".vcf.gz"
        rm_vcfs = ["rm", father_vcf, mother_vcf, child_vcf, new_child_vcf, father_vcf + ".tbi", mother_vcf + ".tbi", child_vcf + ".tbi", new_child_vcf + ".tbi"]
        
        command = bcf_mod + [";"] + bcf_std + [";"] + dng_mod + [";"] + dng_std + [";"] + rm_vcfs
        self.submit_bsub_job(command, job_id, dependent_id=samtools_job_id)
        
        return job_id
    
    def make_bcf_command(self, child, mother, father, bcf, region):
        """ we need to generate BCF for denovogear to work on
        """
        
        region_path = ".{0}.{1}-{2}".format(*region)
        
        # set the paths to the individual VCFs
        child = os.path.splitext(child)[0] + region_path + ".vcf.gz"
        mother = os.path.splitext(mother)[0] + region_path + ".vcf.gz"
        father = os.path.splitext(father)[0] + region_path + ".vcf.gz"
        
        # set up the merge commands, fix the PL field, and generate a BCF for 
        # denovogear
        merge = [self.bcftools, "merge", father, mother, child, "|"]
        pl_fix = ["python", self.pl_fixer, "|"]
        bcf_convert = ["bcftools", "view", "-D", self.dic_path, "-Sb", "-", ">", bcf, ";"]
        
        command = merge + pl_fix + bcf_convert
        
        return command
    
    def construct_dng_command(self, bcf, dnm, region):
        """ Use DNG to call de novo mutations. In this particular case, I use 
        v0.5. As always with DNG, the appropriate PED file (family.ped) has to 
        be in the working directory.
        """
        
        if region[0] != "X":
            dng_chr_type = "auto"
        elif region[0] == "X" and self.proband_sex in self.male_codes:
            dng_chr_type = "XS"
        elif region[0] == "X" and self.proband_sex in self.female_codes:
            dng_chr_type = "XD"
        
        command = [DENOVOGEAR, "dnm", dng_chr_type, "--ped", self.ped_path, "--bcf", bcf, ">", dnm]
        
        return command

def main():
    
    families = open_families(PROBANDS_FILE, FAMILIES_PED_FILE)
    # extract_bams(families, TEMP_DIR)
    
    temp = families[0]
    family = temp[0]
    sex = temp[1]
    
    caller = MosaicCalling(family, sex, TEMP_DIR)
    
    chrom = "1"
    start = "1"
    stop = "5000000"
    region = (chrom, start, stop)
    
    # caller.call_mosaic_de_novos_in_region(region)
    caller.call_mosaic_de_novos()
    
    # for family, sex in families:
        
    #     caller = MosaicCalling(family, sex, TEMP_DIR)
    #     caller.call_mosaic_de_novos()
    


if __name__ == '__main__':
    main()

