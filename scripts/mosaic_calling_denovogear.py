'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from __future__ import division

import os
import argparse
import logging

from mosaic_de_novos.calling import MosaicCalling

logging.basicConfig(filename='mosaic_calling.log',level=logging.DEBUG)

def get_options():
    ''' parse the command line options for the script
    '''
    
    parser = argparse.ArgumentParser(description='Run mosaic denovo calling.')
    parser.add_argument('--proband-bam', required=True, help='BAM file for proband')
    parser.add_argument('--mother-bam', required=True, help='BAM file for mother')
    parser.add_argument('--father-bam', required=True, help='BAM File for father')
    parser.add_argument('--proband-sex', required=True,
        choices=['1', 'M', 'm', 'Male', 'male', '2', 'F', 'f', 'Female', 'female'],
        help='Sex of proband')
    parser.add_argument('--outdir', help='Folder to place denovogear results into')
    
    # and define the region of the genome to call
    parser.add_argument('--chrom', required=True,
        help='Chromosome to find de novos in')
    parser.add_argument('--start',
        help='Region of chromosome to start examining for de novos, omit to '
            'process full chromosome')
    parser.add_argument('--stop',
        help='Region of chromosome to stop examining for de novos, omit to '
            'process full chromosome')
    
    parser.add_argument('--proportion', type=float, default=0.25,
        help='Expected proportion of somatic mosaicism')
    
    parser.add_argument('--generate-merged-bcf', default=False, action='store_true',
        help='Whether to only create a merged BCF for the trio, rather than '
            'running denovogear on a temporary BCF.')
    
    args = parser.parse_args()
    
    if args.proportion > 1 or args.proportion < 0:
        sys.exit('error: argument --proportion: the expected proportion of'
            'somatic mosaicism must be between 0 and 1.')
    
    return args

def main():
    ''' runs mosaic calling for a single region of the genome in a single trio
    '''
    
    args = get_options()
    
    caller = MosaicCalling(args.proband_bam, args.mother_bam, args.father_bam,
        args.proband_sex, args.outdir, args.proportion, args.generate_merged_bcf)
    
    region = (args.chrom, args.start, args.stop)
    
    try:
        caller.call_mosaic_de_novos_in_region(region)
    except KeyboardInterrupt:
        sys.exit(1)
    

if __name__ == '__main__':
    main()
