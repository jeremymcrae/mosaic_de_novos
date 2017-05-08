

import os

BASE = "/lustre/scratch114/projects/ddd/release-main"
BAM_DIRS = [
    "{0}/20140908/sample_improved_bams_hgi_2/".format(BASE), \
    "{0}_2000/20140905/sample_improved_bams_hgi_2/".format(BASE), \
    "{0}_2000/20140910/sample_improved_bams_hgi_2/".format(BASE), \
    "{0}_fy3/20140806/sample_improved_bams_hgi_2/".format(BASE), \
    "{0}_fy3/20140916/sample_improved_bams_hgi_2/".format(BASE), \
    "{0}_y2/20140728/sample_improved_bams_hgi_2/".format(BASE)]


def find_bam_on_lustre(sample_id, SANGER_IDS, all=False):
    """ finds the most appropriate bam to use for a sample.
    
    This depends on the BAMs for samples being available in a few folders.
    
    Args:
        sample_id: DDD sample ID eg DDDP000001
        SANGER_IDS: dictionary of sanger IDs for samples, keyed by DDD stable ID
        all: True/False, whether to return all BAM paths for a sample
    
    Returns:
        path to bam file for sample.
    """
    
    sanger_ids = SANGER_IDS[sample_id]
    
    bam_paths = []
    for folder in BAM_DIRS:
        for sanger_id in sanger_ids:
            path = os.path.join(folder, "{0}.bam".format(sanger_id))
            if os.path.exists(path):
                bam_paths.append(path)
    
    if bam_paths == []:
        return None
    
    if not all:
        # sort the bame paths by creation date, use the most recently created
        bam_paths = sorted(bam_paths, key=os.path.getctime)
        return bam_paths[-1]
    else:
        return bam_paths
