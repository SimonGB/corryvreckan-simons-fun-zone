# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

if [ ! -d "/cvmfs/clicdp.cern.ch" ]; then
    echo "CVMFS not available"
    return
fi

# Get our directory and load the CI init
ABSOLUTE_PATH=`dirname $(readlink -f ${BASH_SOURCE[0]})`

# Set the compiler type to LLVM to also load clang-format and clang-tidy
if [ "$( cat /etc/*-release | grep "Red Hat Enterprise Linux 9" )" ] || [ "$( cat /etc/*-release | grep "AlmaLinux 9" )" ]; then
    export COMPILER_TYPE="llvm"
else
    export COMPILER_TYPE="gcc"
fi


# Load default configuration
source $ABSOLUTE_PATH/../.ci/init_x86_64.sh

# Check if corry executable exists
if [ -f "$ABSOLUTE_PATH/../bin/corry" ]; then
    # Add <path-to-corryvreckan> to PATH
    CORRY_PATH=$( cd "$ABSOLUTE_PATH/../bin" ; pwd -P )
    export PATH=$CORRY_PATH:$PATH
else
    echo "Could not find corry executable. Please complete the installation by executing:"
    echo "$ mkdir build && cd build"
    echo "$ cmake3 .."
    echo "$ make install -j <number_of_cores>"
fi
