""" collects
"""

import os
import sys
import glob
import argparse

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description="Fix PL fields in a VCF file.")
    parser.add_argument("--remove-files", action="store_true", default=False, \
        help="whether to delete the individual denovogear files as we go")
    parser.add_argument("--folder", help="folder containing denovogear output")
    parser.add_argument("--pattern", help="text used to select denovogear files\
        to read from, must be contained within the full path")
    
    args = parser.parse_args()
    
    return args.folder, args.pattern, args.remove_files

def get_denovogear_vars_from_file(path):
    """ extracts variant lines from a denovogear output file
    """
    
    with open(path , "r") as f:
        for line in f:
            if line.startswith("DENOVO"):
                sys.stdout.write(line)

def find_denovogear_paths(folder, path_pattern):
    """ selects denovogear output files by searching for the string pattern
    """
    
    paths = glob.glob(os.path.join(folder, "*denovogear*" + path_pattern + "*.dnm"))
    
    # sort the paths correctly?
    
    return paths

def aggregate_denovogear(folder, path_pattern, remove_files=False):
    """ extract the variants from the individual de novogear files
    """
    
    paths = find_denovogear_paths(folder, path_pattern)
    
    for path in paths:
        get_denovogear_vars_from_file(path)
        if remove_files:
            os.remove(path)

def main():
    
    path, path_pattern, remove_files = get_options()
    de_novos = aggregate_denovogear(path, path_pattern, remove_files)

if __name__ == '__main__':
    main()

