# SPDX-FileCopyrightText: 2022-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

import csv
import os.path
import logging
import copy

import sys
from sys import exit # use sys.exit instead of built-in exit (latter raises exception)

import misc

class CommentedFile:
    """ Decorator for text files: filters out comments (i.e. first char of line #)
    Based on http://www.mfasold.net/blog/2010/02/python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines/

    """
    def __init__(self, f, comment_string="#"):
        self.f = f
        self.comment_string = comment_string
        self.line_count = 0
    def rewind(self):
        self.f.seek(0)
        self.line_count = 0
    def __next__(self):
        line = self.f.__next__()
        self.line_count += 1
        while line.startswith(self.comment_string) or not line.strip(): # test if line commented or empty
            line = self.f.__next__()
            self.line_count += 1
        return str(line)
    def __iter__(self):
        return self

class Loader():
    """ Load and parse the csv file for the given set of runs and
    return nested dictionary: a collection of dictionaries, one for
    each csv row matching a run number.

    """

    def __init__(self, filename, runs):
        self.filename = filename
        self.runs = runs
        self.log = logging.getLogger('jobsub')

        self.parameters = {} # store all information needed from the csv file

    def parse(self):
        if self.filename is None:
            return False # if no file name given, return empty collection here
        if not os.path.isfile(self.filename): # check if file exists
            self.log.error("Could not find the specified csv file '"+self.filename+"'!")
            exit(1)
        try:
            self.log.debug("Opening csv file '"+self.filename+"'.")
            csv_file = open(self.filename, 'rt')
            filtered_file = CommentedFile(csv_file)

            try:
                # construct a sample for the csv format sniffer:
                sample = ""
                try:
                    while (len(sample)<1024):
                        sample += filtered_file.__next__()
                except StopIteration:
                    self.log.debug("End of csv file reached, sample limited to " + str(len(sample))+ " bytes")
                dialect = csv.Sniffer().sniff(sample) # test csv file format details
                dialect.escapechar = "\\"
                self.log.debug("Determined the CSV dialect as follows: delimiter=%s, doublequote=%s, escapechar=%s, lineterminator=%s, quotechar=%s , quoting=%s, skipinitialspace=%s", dialect.delimiter, dialect.doublequote, dialect.escapechar, list(ord(c) for c in dialect.lineterminator), dialect.quotechar, dialect.quoting, dialect.skipinitialspace)
                filtered_file.rewind() # back to beginning of file
                reader = csv.DictReader(filtered_file, dialect=dialect) # now process CSV file contents here and load them into memory
                reader.__next__() # python requires an actual read access before filling 'DictReader.fieldnames'

                self.log.debug("CSV file contains the header info: %s", reader.fieldnames)
                try:
                    reader.fieldnames = [field.lower() for field in reader.fieldnames] # convert to lower case keys to avoid confusion
                    reader.fieldnames = [field.strip() for field in reader.fieldnames] # remove leading and trailing white space
                except TypeError:
                    self.log.error("Could not process the CSV file header information. csv.DictReader returned fieldnames: %s", reader.fieldnames)
                    exit(1)
                if not "runnumber" in reader.fieldnames: # verify that we have a column "runnumber"
                    self.log.error("Could not find a column with header label 'RunNumber' in file '"+self.filename+"'!")
                    exit(1)
                if "" in reader.fieldnames:
                    self.log.warning("Column without header label encountered in csv file '"+self.filename+"'!")
                self.log.info("Successfully loaded csv file'"+self.filename+"'.")
                # first: search through csv file to find corresponding runnumber entry line for every run
                filtered_file.rewind() # back to beginning of file
                reader.__next__()   # .. and skip the header line
                missing_runs = list(self.runs) # list of runs to look for in csv file

                cnt = 0
                for row in reader: # loop over all rows once
                    try:
                        cnt += 1
                        runs = [int(val) for val in misc.parseBrackets(row["runnumber"])]
                        missing_runs = [run for run in missing_runs if not run in runs]
                        # self.log.debug(f"Found entry in csv file for run {run} on line {filtered_file.line_count}")

                        self.parameters[cnt-1] = {} # start counting at 0
                        self.parameters[cnt-1].update(row) # start counting at 0
                    except ValueError: # int conversion error
                        self.log.warn("Could not interpret run number on line "+str(filtered_file.line_count)+" in file '"+self.filename+"'.")
                        continue
                if len(missing_runs)==0:
                    self.log.debug("Found at least one line for each run we were searching for.")

                self.log.debug("Searched over "+str(filtered_file.line_count)+" lines in file '"+self.filename+"'.")
                if not len(missing_runs)==0:
                    self.log.error("Could not find an entry for the following run numbers in '"+self.filename+"': "+', '.join(map(str, missing_runs)))
            finally:
                csv_file.close()
        except csv.Error as e:
            self.log.error("Problem loading the csv file '"+self.filename+"': %s"%e)
            exit(1)

        return True

    def process(self, runnr, steering_string_base):
        appendixes = []
        steering_strings = []
        for line in self.parameters: # go through line by line
            # if we have a csv file we can parse, we will check for the runnumber and replace any
            # variables identified by the csv header by the run specific value
            try:
                # check if the runnumber in csv is provided as list of values or range
                runnrs = misc.parseBrackets(self.parameters[line]["runnumber"])
                if not int(runnr) in runnrs and not runnr in runnrs:
                    continue
                # process only the selected runnr
                parameters = copy.copy(self.parameters[line])
                parameters["runnumber"] = runnr

                results = []
                self.log.debug(f"Assembling parameters for run: {runnr}")
                self.parameterAssembler(parameters, results)

                for parameter_combination in results:
                    self.log.debug(f"Parameter combination: {parameter_combination}")
                    # make a copy of the preprocessed steering file content
                    steering_strings.append(steering_string_base)
                    appendixes.append('_' + str(parameter_combination["runnumber"]))

                    for key in parameter_combination.keys():
                        # check if we actually find all parameters from the csv file in the steering file - warn if not
                        self.log.debug(f"Parsing steering file for csv field name {key}")
                        try:
                            # check that the field name is not empty and do not yet replace the runnumber
                            if not key == "" and not key == "runnumber":
                                steering_strings[-1] = misc.ireplace("@" + key + "@", parameter_combination[key], steering_strings[-1])
                                if parameter_combination[key] == "@RunNumber@":
                                    appendix = str(parameter_combination["runnumber"])
                                else:
                                    appendix = str(parameter_combination[key])
                                appendixes[-1] = appendixes[-1] + '-' + key + ':' + appendix
                                self.log.debug("appendix is now '%s'", appendixes[-1])
                        except EOFError:
                            self.log.warn("Parameter '" + key + "' from the csv file was not found in the template file (already overwritten by config file parameters?)")
                    # finally replace runnumber - this is in case one of the parameters in options is matching runnumber
                    steering_strings[-1] = misc.ireplace("@RunNumber@", parameter_combination["runnumber"], steering_strings[-1])
            except KeyError:
                self.log.warning("Run #" + runnr + " was not found in the specified CSV file - will skip this run! ")
                continue

        return appendixes, steering_strings

    def parameterAssembler(self, parameters, result, selection = {}, key_iterator = None):
        current_parameter = []
        if key_iterator is None:
            key_iterator = iter(parameters)
        try:
            key = next(key_iterator)
            self.log.debug(f"Parsing key: {key}")
        except StopIteration:
            raise StopIteration(f'Finished dictionary {parameters}')

        # remove all whitespaces from beginning and end of string (not in the middle)
        parameters[key] = parameters[key].strip()

        # process list or range
        current_parameter = misc.parseBrackets(parameters[key])

        # recursion to access all possible combinations
        for parameter in current_parameter:
            selection[key] = parameter
            try:
                key_iterator_next = copy.copy(key_iterator)
                self.parameterAssembler(parameters, result, selection, key_iterator_next)
            except StopIteration:
                self.log.debug(f"Reached last key: {key}, appending parameter selection: {selection}")
                result.append(copy.copy(selection))
