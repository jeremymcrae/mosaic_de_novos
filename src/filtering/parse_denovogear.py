""" class to open denovogear output
"""

class ParseDenovogear(object):
    """ class to read denovogear output
    """
    
    def __init__(self, path):
        """ initiates the class with the path and reads the file
        
        Args:
            path: path to a denovogear output file
        """
        
        self.variants = {}
        self._open_file(path)
     
    def _open_file(self, path):
        """ opens denovogear output, and parses variant lines
        
        Args:
            path: path to a denovogear output file
        """
        
        with open(path, "r") as f:
            for line in f:
                if line.startswith("DENOVO"):
                    var = self._parse_line(line)
                    chrom = self._convert_chrom_to_numeric(var["ref_name"])
                    key = (chrom, int(var["coor"]))
                    self.variants[key] = var
    
    def _convert_chrom_to_numeric(self, chrom):
        """ converts a chromosome string to equivalent integer, for easy sorting
        
        Args:
            chrom: string of the chromosome ("ref_name" in denovogear output)
        """
        
        sex_chroms = {"X": 23, "Y": 24}
        
        try:
            chrom = int(chrom)
        except ValueError:
            if chrom in sex_chroms:
                chrom = sex_chroms[chrom]
        
        return chrom
    
    def _parse_line(self, line):
        """ parses a denovogear variant output line
        
        Args:
            line: line for a variant in a denovogear output fields
        
        Returns:
            a dictionary filled with entries from the denovogear line
        """
        
        fields = line.strip().split(" ")
        read_idx = fields.index("READ_DEPTH")
        map_idx = fields.index("MAPPING_QUALITY")
        
        variant = {}
        variant["type"] = fields[0]
        
        # grab all the fields prior to the read depth fields, since these are
        # easy unique key, value pairs up until then.
        for i in range(2, read_idx, 2):
            key = fields[i].strip(":")
            variant[key] = fields[i + 1]
        
        # and parse the read depth and quality metrics
        variant = self._parse_read_depth(fields, read_idx, variant)
        variant = self._parse_qual(fields, read_idx, variant)
        
        # might as well hang on to the original line, if we want to simply dump
        # these for later analysis
        variant["line"] = line
        
        # a very few lines have the ID as "id", not "ID". Standardise this field
        if "id" in variant:
            variant["ID"] = variant["id"]
            del variant["id"]
        
        variant["alt_allele"] = self._parse_alt_allele(variant)
        if len(variant["alt_allele"]) == 1:
            variant["alt_allele"] = variant["alt_allele"].pop()
        
        return variant
    
    def _parse_read_depth(self, split_line, index, variant):
        """ grab the person-specific read depth fields,
        
        The read depth fields are unfortunately named the same as for the qual
        fields, so we need to read the line, starting from an index position.
        
        Args:
            split_line: denovogear variant line
            index: position in the split line where the read depth fields start
            variant: dictionary entry for the current variant
        
        Returns:
            variant, with read depth entries added
        """
        
        variant["READ_DEPTH"] = {}
        for i in range(index + 1, index + 6, 2):
            key = split_line[i].strip(":")
            variant["READ_DEPTH"][key] = split_line[i + 1]
        
        return variant
        
    def _parse_qual(self, split_line, index, variant):
        """ grab the person-specific mapping quality fields,
        
        The qual fields are unfortunately named the same as for the read depth
        fields, so we need to read the line, starting from an index position.
        
        Args:
            split_line: denovogear variant line
            index: position in the split line where the qual fields start
            variant: dictionary entry for the current variant
        
        Returns:
            variant, with qual entries added
        """
        
        variant["MAPPING_QUALITY"] = {}
        for i in range(index + 1, index + 6, 2):
            key = split_line[i].strip(":")
            variant["MAPPING_QUALITY"][key] = split_line[i + 1]
        
        return variant
    
    def _parse_alt_allele(self, variant):
        """ grab the alt allele for the variant
        """
        
        alt_allele = set(variant["ALT"].split(","))
        alt_allele = alt_allele - set(["X"])
        
        # get the alleles from the targetted genotype
        alleles = set(list(variant["tgt"].split("/")[0]))
        
        # intersect the alt alleles with the alleles in the targetted genotypes
        if len(alt_allele) > 1:
            alt_allele = alt_allele & alleles
        
        # make sure we only have a single alt allele left
        alt_allele = list(alt_allele)
        
        return alt_allele
     
    def __sub__(self, other):
        """ get the variants that are unique to the current instance
        """
        
        return set(self.variants) - set(other.variants)
    
    def __iter__(self):
        """ get an iterable for the class, to iterate through the variants
        """
        
        return iter(self.variants)
    
    def __getitem__(self, key):
        """ gets a variant from the class
        
        Args:
            key: (chrom, pos) tuple,  which must occur in self.variants
        
        Returns:
            dictionary entry for variant with matching key
        """
        
        return self.variants[key]
