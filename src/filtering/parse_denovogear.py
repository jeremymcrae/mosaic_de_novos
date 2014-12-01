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
                    key = (var["ref_name"], var["coor"])
                    self.variants[key] = var
    
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
        
        # and parse the 
        variant = self._parse_read_depth(fields, read_idx, variant)
        variant = self._parse_qual(fields, read_idx, variant)
        
        # might as well hang on to the original line, if we want to simply dump 
        # these for later analysis
        variant["line"] = line
        
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
    
    def get_variants(self):
        """ get the variants loaded from a denovogear file
        
        Returns:
            dictionary of parsed variant lines, indexed by (chrom, pos) tuples
        """
        
        return self.variants
     
    def __sub__(self, other):
        """ get the variants that are unique to the current instance
        """
        
        return set(self.get_variants()) - set(other.get_variants())
    
    def __iter__(self):
        return iter(self.variants)
    
    def get_subset(self, keys):
        """ gets a subset of variants given a list of keys
        
        Args:
            keys: list of (chrom, pos) tuples, all of which must occur in
                self.variants
        
        Returns:
            dictionary of variants, restricted to the set with matching keys
        """
        
        variants = self.get_variants()
        
        subset = {}
        for key in keys:
            if key in variants:
                subset[key] = variants[key]
        
        return subset


