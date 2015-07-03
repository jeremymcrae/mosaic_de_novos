""" functions to submit compute jobs to a LSF cluster, and check for submitted
jobs.
"""

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None, queue="normal", cpus=1):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    threads=""
    if cpus >1:
        threads="-n{0} -R 'span[hosts=1]'".format(cpus)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", queue, "-o", log, mem, threads]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string(prefix=None):
    """ make a random string, which we can use for bsub job IDs, so that
    different jobs do not have the same job IDs.
    """
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    # done't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    if prefix is not None:
        hash_string = prefix + hash_string
    
    return hash_string

def get_bjobs():
    """ get a list of submitted jobs
    """
    
    command = ["bjobs", "-o", "\"JOBID", "USER", "STAT", "QUEUE", "JOB_NAME", "delimiter=';'\""]
    command = " ".join(command)
    output = subprocess.check_output(command, shell=True, stderr=open(os.devnull, "w"))
    
    bjobs = []
    for line in output.split("\n"):
        if line.startswith("JOBID") or line == "":
            continue
        
        line = line.strip().split(";")
        entry = {"jobid": line[0], "user":line[1], "stat":line[2], \
            "queue":line[3], "job_name":line[4]}
        bjobs.append(entry)
    
    return bjobs
