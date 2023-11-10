# SPDX-FileCopyrightText: 2022-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

import logging

def parseIntegerString(inputstr=""):
    """
    return a list of selected values when a string in the form:
    1-4,6
    would return:
    1,2,3,4,6
    as expected...
    (from http://thoughtsbyclayg.blogspot.de/2008/10/parsing-list-of-numbers-in-python.html)

    Modified such that it returns a list of strings
    if the conversion to integer fails, e.g.
    "10ns, 20ns"
    would return:
    "10ns", "20ns"
    """
    selection = list()
    # tokens are comma separated values
    tokens = [substring.strip() for substring in inputstr.split(',')]
    for i in tokens:
        try:
            # typically tokens are plain old integers
            selection.append(int(i))
        except ValueError:
            try:
                # if not, then it might be a range
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items separated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token)-1]
                    for value in range(first, last+1):
                        selection.append(value)
            except ValueError:
                # if not treat as string, not integer
                selection.append(i)
    return selection # end parseIntegerString

def parseBrackets(inputstr):
    log = logging.getLogger('jobsub')
    selection = []
    if inputstr[0] == '{':
        log.debug("Found open bracket, look for matching close bracket.")
        if inputstr[-1] == '}':
            log.debug("Found matching close bracket, Interpret as range or set of parameters.")
            # remove curly brackets:
            substring = inputstr.strip("{}")
            # Check if string contains "," or "-", i.e. a set or range of values
            # If not, no conversion is required (or even possible in case of file paths etc.)
            # If yes, call parseIntegerString() and list with all values.
            if any(delimiter in substring for delimiter in [',','-']):
                selection = parseIntegerString(substring)
                log.debug(f"Found delimiter for {substring}, produced {selection}")
            else:
                # current_parameter needs to be a list to get len(list) = 1
                selection.append(inputstr)
                log.debug(f"No delimiter found for {substring}")
        else:
            log.error(f"No matching close bracket found. Please update the string {inputstr}.")
            exit(1)
    else:
        log.debug(f"No bracket found, interpret as one string {inputstr}.")
        selection.append(inputstr)
    return selection


def ireplace(old, new, text):
    """
    case insensitive search and replace function searching through string and returning the filtered string
    (based on http://stackoverflow.com/a/4773614)

    """
    idx = 0
    occur = 0
    while idx < len(text):
        index_l = text.lower().find(old.lower(), idx)
        if index_l == -1:
            if occur == 0:
                raise EOFError("Could not find string "+old)
            return text
        text = text[:index_l] + str(new) + text[index_l + len(old):]
        idx = index_l + len(str(new))
        occur = occur+1
    if occur == 0:
        raise EOFError("Could not find string "+old)
    return text

def checkProgram(name):
    """ Searches PATH environment variable for executable given by parameter """
    import os
    for dir in os.environ['PATH'].split(os.pathsep):
        prog = os.path.join(dir, name)
        if os.path.exists(prog): return prog

def checkSteer(sstring):
    """ Check string for any occurrence of @.*@ and return boolean. """
    log = logging.getLogger('jobsub')
    import re
    hits = re.findall("@.*@", sstring)
    if hits:
        log.error ("Missing configuration parameters: "+', '.join(map(str, hits)))
        return False
    else:
        return True

def createSteeringFile(log, args, steering_string, suffix):
    """ Create file with replaceable parameters defined by @.*@ and its name. """
    import os

    if not checkSteer(steering_string):
        raise ValueError("Missing steeringString " + steering_string)

    # update this line too
    log.debug ("Writing steering file for run " + suffix)

    # Get "jobtask" as basename of the configuration file:
    jobtask = os.path.splitext(os.path.basename(args.conf_file))[0]
    # Write the steering file:
    filename = jobtask + "_" + suffix
    log.info("filename = " + filename)
    steering_file = open(filename+".conf", "w")

    try:
        steering_file.write(steering_string)
    finally:
        steering_file.close()

    return filename

def zipLogs(path, filename):
    """  stores output from Corryvreckan in zip file; enables compression if necessary module is available """
    import zipfile
    import os.path
    log = logging.getLogger('jobsub')
    try:     # compression module might not be available, therefore try import here
        import zlib
        compression = zipfile.ZIP_DEFLATED
        log.debug("Creating *compressed* log archive")
    except ImportError: # no compression module available, use flat files
        compression = zipfile.ZIP_STORED
        log.debug("Creating flat log archive")
    try:
        zf = zipfile.ZipFile(os.path.join(path, filename)+".zip", mode='w') # create new zip file
        try:
            zf.write(os.path.join("./", filename)+".conf", compress_type=compression) # store in zip file
            zf.write(os.path.join("./", filename)+".log", compress_type=compression) # store in zip file
            os.remove(os.path.join("./", filename)+".conf") # delete file
            os.remove(os.path.join("./", filename)+".log") # delete file
            log.info("Logs written to "+os.path.join(path, filename)+".zip")
        finally:
            log.debug("Closing log archive file")
            zf.close()
    except IOError: # could not create zip file - path non-existent?!
        log.error("Input/Output error: Could not create log and steering file archive ("+os.path.join(path, filename)+".zip"+")!")

def poolChecker(results, heartbeat = 1):
    """  checks the status of parallel pool """
    import time
    log = logging.getLogger('jobsub')
    if results == []:
        log.warning("There were problems with the submission")
    elif None not in results:
        # parallel loop checker from https://stackoverflow.com/a/70666333
        while True:
            log.debug("Heartbeat")
            time.sleep(heartbeat)
            # catch exception if results are not ready yet
            try:
                ready = [result.ready() for result in results]
                successful = [result.successful() for result in results]
                rcodes = [result.get() for result in results]
            except Exception as e:
                # not_ready exception is always thrown, ignore it
                continue
            # exit loop if all tasks returned success
            if all(successful):
                if any(rcodes):
                    log.warning(f"There was a problem with submission")
                else:
                    log.info(f"All job submissions finished without problems")
                break
            # raise exception reporting exceptions received from workers
            if all(ready) and not all(successful):
                raise Exception(f'Workers raised following exceptions {[result._value for result in results if not result.successful()]}')
