# SPDX-FileCopyrightText: 2022-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

# need some additional libraries for process interaction
import datetime
import asyncio
import os
import sys
from asyncio.subprocess import PIPE

import shlex

# runner class for local job submission
class Runner():
    def __init__(self, log, filename, silent):
        self.log = log
        self.silent = silent
        self.filename = filename

    # Based on the solution proposed here:
    # https://gitlab.cern.ch/corryvreckan/corryvreckan/-/issues/146#note_4465941
    async def read_stream_and_display(self, stream, display):
        """Read from stream line by line until EOF, display, and capture the lines."""
        while True:
            line = await stream.readline()
            if not line:
                break

            if not self.silent:
                if b'WARNING' in line.strip():
                    self.log.warning(line.strip())
                elif b'ERROR' in line.strip():
                    self.log.error(line.strip())
                elif b'FATAL' in line.strip():
                    self.log.critical(line.strip())
                else:
                    self.log.info(line.strip())
            self.log_file.write(str(line,'utf-8'))

    async def read_and_display(self, cmd):
        """Capture cmd's stdout, stderr while displaying them as they arrive (line by line)."""
        # start process
        self.log.info("Starting process %s", cmd)
        process = await asyncio.create_subprocess_exec(*shlex.split(cmd), stdout=PIPE, stderr=PIPE)

        # read child's stdout/stderr concurrently (capture and display)
        try:
            # stdout, stderr = await asyncio.gather(
            await asyncio.gather(
                self.read_stream_and_display(process.stdout, sys.stdout.buffer.write),
                self.read_stream_and_display(process.stderr, sys.stderr.buffer.write))
        except Exception:
            process.kill()
            raise
        finally:
            # wait for the process to exit
            rc = await process.wait()

        self.log.info("Finished process %s", cmd)
        return rc

    def run(self, cmd):
        # open log file
        self.log_file = open(self.filename+".log", "w")
        # print timestamp to log file
        self.log_file.write("---=== Analysis started on " + datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p") + " ===---\n\n")

        # run process
        try:
            # execute the code and monitor outputs
            rcode = asyncio.run(self.read_and_display(cmd))
        except Exception:
            raise
        finally:
            self.log_file.write("---=== Analysis finished on " + datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p") + " ===---\n\n")
            self.log_file.close()

        return rcode
