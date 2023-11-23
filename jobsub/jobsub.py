#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2022-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

"""
jobsub: a tool for Corryvreckan job submission

Steering files are generated based on steering templates by substituting
job and run specific variables. The variable information can be
provided to jobsub by command line argument, by a config file or by a
text file with run parameters (in comma-separated value/csv format).

Run
python jobsub.py --help
to see the list of command line options.

"""
import argparse
import sys
import logging
import misc
import os
import multiprocessing
import hashlib
import signal

def runCorryvreckanLocally(filename, jobtask, silent):
    """ Runs Corryvreckan and stores log of output """
    log = logging.getLogger('jobsub.' + jobtask)

    rcode = None # the return code that will be set by a later subprocess method

    # check for Corryvreckan executable
    cmd = misc.checkProgram("corry")
    if cmd:
        log.debug("Found Corryvreckan executable: " + cmd)
    else:
        log.error("Corryvreckan executable not found in PATH!")
        log.error(os.getcwd())
        return 1

    # search for stdbuf command: adjust stdout buffering
    stdbuf = misc.checkProgram("stdbuf")
    if stdbuf:
        log.debug("Found stdbuf, will use line buffered output.")
        # -oL: adjust standard output stream buffering to line buffered
        cmd = stdbuf + " -oL " + cmd

    import runner

    cmd = cmd+" -c "+filename+".conf"
    try:
        # execute the asynchronous job submission
        r = runner.Runner(log, filename, silent)
        rcode = r.run(cmd)
    except OSError as e:
        log.critical("Problem with Corryvreckan execution: Command '%s' resulted in error %s", cmd, e)
        return 2
    return rcode

def runCorryvreckanCondor(filename, subfile, jobtask):
    """ Submits the Corryvreckan job to HTCondor """
    log = logging.getLogger('jobsub.' + jobtask)
    # We are running on HTCondor.

    rcode = None

    # check for qsub executable
    cmd = misc.checkProgram("condor_submit")
    if cmd:
        log.debug("Found condor_submit executable: " + cmd)
    else:
        log.error("condor_submit executable not found in PATH!")
        return 1

    # Add condor_submit parameters:
    cmd = cmd+" -batch-name \"Corry"+jobtask+"\" "

    # check for Corryvreckan executable
    corry = misc.checkProgram("corry")
    if corry:
        log.debug("Found Corryvreckan executable: " + corry)
        cmd = cmd+" executable="+corry
    else:
        log.error("Corryvreckan executable not found in PATH!")
        return 1

    cmd = cmd+" arguments=\"-c "+os.path.abspath(filename+".conf")+"\""

    # Add Condor submission configuration file:
    cmd = cmd+" "+subfile

    try:
        # run process
        log.info ("Now submitting Corryvreckan job: "+filename+".conf to HTCondor")
        log.debug ("Executing: "+cmd)
        os.popen(cmd)
        rcode = 0
    except OSError as e:
        log.critical("Problem with HTCondor submission: Command '%s' resulted in error %s", cmd, e)
        return 2
    return rcode

def submitJobs(log, pool, args, filename, parameters):
    # bail out if running a dry run
    rcode = None
    if args.dry_run:
        log.info("Dry run: skipping Corryvreckan execution. Steering file written to "+filename+'.conf')
    elif args.htcondor_file:
        rcode = pool.apply_async(runCorryvreckanCondor, (filename, args.htcondor_file, filename))
        #rcode = runCorryvreckanCondor(filename, args.htcondor_file, filename) # start HTCondor submission
        #if rcode == 0:
        log.info("HTCondor job submitted")
        #else:
        #    log.error("HTCondor submission returned with error code "+str(rcode))
    else:
        rcode = pool.apply_async(runCorryvreckanLocally, (filename, filename, args.silent))
        #rcode = runCorryvreckanLocally(filename, filename, args.silent) # start Corryvreckan execution
        #if rcode == 0:
        log.info("Local job submitted")
        #else:
        #    log.error("Corryvreckan returned with error code "+str(rcode))
        #misc.zipLogs(parameters["logpath"], filename)
    return rcode

def main(argv=None):
    """  main routine of jobsub: a tool for Corryvreckan job submission """
    log = logging.getLogger('jobsub') # set up logging
    formatter = logging.Formatter('%(name)s(%(levelname)s): %(message)s',"%H:%M:%S")
    handler_stream = logging.StreamHandler()
    handler_stream.setFormatter(formatter)
    log.addHandler(handler_stream)
    # using this decorator, we can count the number of error messages
    class callcounted(object):
        """Decorator to determine number of calls for a method"""
        def __init__(self,method):
            self.method=method
            self.counter=0
        def __call__(self,*args,**kwargs):
            self.counter+=1
            return self.method(*args,**kwargs)
    log.error=callcounted(log.error)

    if argv is None:
        argv = sys.argv
        progName = os.path.basename(argv.pop(0))

    # command line argument parsing
    parser = argparse.ArgumentParser(prog=progName, description="A tool for the convenient run-specific modification of Corryvreckan configuration files and their execution through the corry executable")
    parser.add_argument('--version', action='version', version='Revision: $Revision$, $LastChangedDate$')
    parser.add_argument("-c", "--conf-file", "--config", help="Configuration file with all Corryvreckan algorithms defined", metavar="FILE")
    parser.add_argument('--option', '-o', action='append', metavar="NAME=VALUE", help="Specify further options such as 'beamenergy=5.3'. This switch be specified several times for multiple options or can parse a comma-separated list of options. This switch overrides any config file options and also overwrites hard-coded settings on the Corryvreckan configuration file.")
    parser.add_argument("-htc", "--htcondor-file", "--batch", help="Specify condor_submit parameter file for HTCondor submission. Run HTCondor submission via condor_submit instead of calling Corryvreckan directly", metavar="FILE")
    parser.add_argument("-csv", "--csv-file", help="Load additional run-specific variables from table (text file in csv format)", metavar="FILE")
    parser.add_argument("--log-file", help="Save submission log to specified file", metavar="FILE")
    parser.add_argument("-v", "--verbosity", default="info", help="Sets the verbosity of log messages during job submission where LEVEL is either debug, info, warning or error", metavar="LEVEL")
    parser.add_argument("-s", "--silent", action="store_true", default=False, help="Suppress non-error (stdout) Corryvreckan output to console")
    parser.add_argument("--dry-run", action="store_true", default=False, help="Write configuration files but skip actual Corryvreckan execution")
    parser.add_argument("--subdir", action="store_true", default=False, help="Execute every job in its own subdirectory instead of all in the base path")
    parser.add_argument("--plain", action="store_true", default=False, help="Output written to stdout/stderr and log file in prefix-less format i.e. without time stamping")
    parser.add_argument("-j", "--cores", metavar='N', type=int, default=1, help="Number of cores used for the local job submission")
    parser.add_argument("--zfill", metavar='N', type=int, help="Fill run number with zeros up to the defined number of digits")
    parser.add_argument("runs", help="The runs to be analyzed; can be a list of single runs and/or a range, e.g. 1056-1060.", nargs='*')
    args = parser.parse_args(argv)

    # Try to import the colorer module
    try:
        import Colorer
    except ImportError:
        pass

    # set the logging level
    numeric_level = getattr(logging, "INFO", None) # default: INFO messages and above
    if args.verbosity:
        # Convert log level to upper case to allow the user to specify --log=DEBUG or --log=debug
        numeric_level = getattr(logging, args.verbosity.upper(), None)
        if not isinstance(numeric_level, int):
            log.error('Invalid log level: %s' % args.verbosity)
            return 2
    handler_stream.setLevel(numeric_level)
    log.setLevel(numeric_level)

    if args.plain:
        formatter = logging.Formatter('%(message)s')
        handler_stream.setFormatter(formatter)

    # set up submission log file if requested on command line
    if args.log_file:
        handler_file = logging.FileHandler([args.log_file])
        handler_file.setFormatter(formatter)
        handler_file.setLevel(numeric_level)
        log.addHandler(handler_file)

    # check existence of htcondor file
    if args.htcondor_file:
        args.htcondor_file = os.path.abspath(args.htcondor_file)
        if not os.path.isfile(args.htcondor_file):
            log.critical("HTCondor submission parameters file '"+args.htcondor_file+"' not found!")
            return 1

    log.debug( "Command line arguments used: %s ", args )

    runs = list()
    for runnum in args.runs:
        try:
            log.debug("Parsing run-range argument: '%s'", runnum)
            runs = runs + misc.parseIntegerString(runnum)
        except ValueError:
            log.error("The list of runs contains non-integer and non-range values: '%s'", runnum)
            return 2

    if not runs:
        log.error("No run numbers were specified. Please see '"+progName+" --help' for details.")
        return 2

    if len(runs) > len(set(runs)): # sets items are unique
        log.error("At least one run is specified multiple times!")
        return 2

    # dictionary keeping our parameters
    # here you can set some minimal default config values that will (possibly) be overwritten by the config file
    parameters = {"logpath":"."}

    # Parse option part of the  argument here -> overwriting config options
    if args.option is None:
        log.debug("Nothing to parse: No additional config options specified through command line arguments. ")
    else:
        try:
            # now parse any options given through the -o cmd line switch
            cmdoptions = dict(opt.strip().split('=', 1) for optlist in args.option for opt in optlist.split(',')) # args.option is a list of lists of strings we need to split at every '='
        except ValueError:
            log.error( "Command line error: cannot parse --option argument(s). Please use a '--option name=value' format. ")
            return 2
        for key in cmdoptions: # and overwrite our current config settings
            log.debug( "Parsing cmd line: Setting "+key+" to value '"+cmdoptions[key]+"', possibly overwriting corresponding config file option")
            parameters[key.lower()] = cmdoptions[key]

    log.debug( "Our final config:")
    for key, value in parameters.items():
        log.debug ( "     "+key+" = "+value)

    if not os.path.isfile(args.conf_file):
        log.critical("Configuration template '"+args.conf_file+"' not found!")
        return 1

    log.debug( "Opening configuration template "+args.conf_file)
    steering_string_base = open(args.conf_file, "r").read()

    #Query replace steering template with our parameter set
    log.debug ("Generating base configuration file")
    for key in parameters.keys():
        # check if we actually find all parameters from the config in the steering file
        try:
            steering_string_base = misc.ireplace("@" + key + "@", parameters[key], steering_string_base)
        except EOFError:
            if (not key == "conf_file" and not key == "logpath"): # do not warn about default content of config
                log.warning("Parameter '" + key + "' was not found in configuration template "+args.conf_file)

    # CSV table
    if args.csv_file:
        log.debug ("Loading csv file (if requested)")
        import loader
        csv_loader = loader.Loader(args.csv_file, runs)
        csv_parsed = csv_loader.parse()

    # Set signal to ignore SIGINT before spawning the pool workers
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    # execute all submissions in a parallel fashion
    pool = multiprocessing.Pool(processes=args.cores)

    # setup mechanism to deal with user pressing ctrl-c in a safe way while we execute Corryvreckan later
    keep_running = {'Sigint':'no'}
    def signal_handler(signal, frame):
        """ log if SIGINT detected, set variable to indicate status """
        log.critical ('You pressed Ctrl+C!')
        pool.close() # prevent any further additions to the pool
        keep_running['Sigint'] = 'seen'
    prevINTHandler = signal.signal(signal.SIGINT, signal_handler)

    # ----------------------------------------------------------------
    # Actual processing and submission of the runs
    log.info("Will now start processing the following runs: "+', '.join(map(str, runs)))
    results = []
    for run in runs:
        if keep_running['Sigint'] == 'seen':
            log.critical("Stopping to process remaining runs now")
            break  # if we received ctrl-c (SIGINT) we stop processing here

        if args.zfill:
            runnr = str(run).zfill(args.zfill)
        else:
            runnr = str(run)
        log.info ("Now generating configuration file for run number "+runnr+"..")

        # When  running in subdirectories for every job, create it:
        if args.subdir:
            run_path = "run_"+runnr
            if not os.path.exists(run_path):
                os.makedirs(run_path)

            # Descend into subdirectory:
            base_path = os.getcwd()
            os.chdir(run_path)

        suffixes = []
        steering_strings = []
        if args.csv_file and csv_parsed:
            # runs only if csv properly parsed
            suffixes, steering_strings = csv_loader.process(runnr, steering_string_base)
        else:
            # Used for CLI parsing
            try:
                suffixes.append(runnr)
                steering_strings.append(misc.ireplace("@RunNumber@", str(runnr), steering_string_base))
            except KeyError:
                log.error("Could not change RunNumber in the config file")
                continue

        # ----------------------------------------------------------------
        # Job submission based on selected method - local/condor
        for suffix, steering_string in zip(suffixes, steering_strings):
            try:
                suffix_hash = f'{run}_' + hashlib.shake_128(suffix.encode()).hexdigest(16)
                steering_filename = misc.createSteeringFile(log, args, steering_string, suffix_hash)
                results.append(submitJobs(log, pool, args, steering_filename, parameters))
            except Exception as e:
                log.error(f"Could not create submission file with suffix {suffix_hash} due to {e}")

        # Return to old directory:
        if args.subdir:
            os.chdir(base_path)

    misc.poolChecker(results)

    signal.signal(signal.SIGINT, prevINTHandler)
    if log.error.counter>0:
        log.warning("There were "+str(log.error.counter)+" error messages reported")

    return 0

if __name__ == "__main__":
    sys.exit(main())
