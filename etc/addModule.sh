#!/bin/bash

# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

echo -e "\nPreparing code basis for a new module:\n"

# Ask for module name:
read -p "Name of the module? " MODNAME

# Ask for module type:
echo -e "Type of the module?\n"
type=0
select yn in "global" "detector" "dut"; do
    case $yn in
        global ) type=1; break;;
        detector ) type=2; break;;
        dut ) type=3; break;;
    esac
done

# If type is detector, ask if aux should be excluded:
exclude_aux=0
if [ "$type" == 2 ]; then
    echo -e "Should auxiliary detectors be excluded from this module?\n"
    select yn in "yes" "nope"; do
        case $yn in
            yes ) exclude_aux=1; break;;
            nope ) exclude_aux=0; break;;
        esac
    done
fi

# If type is detector, ask if pass should be excluded:
exclude_pass=1
if [ "$type" == 2 ]; then
    echo -e "Should passive detectors be excluded from this module?\n"
    select yn in "yes" "nope"; do
        case $yn in
            yes ) exclude_pass=1; break;;
            nope ) exclude_pass=0; break;;
        esac
    done
fi

echo "Creating directory and files..."

echo
# Try to find the modules directory:
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
DIRECTORIES[0]="${BASEDIR}/../../src/modules"
DIRECTORIES[1]="../src/modules"
DIRECTORIES[2]="src/modules"

MODDIR=""
for DIR in "${DIRECTORIES[@]}"; do
    if [ -d "$DIR" ]; then
        MODDIR="$DIR"
        break
    fi
done

# Create directory
mkdir "$MODDIR/$MODNAME"

# Copy over CMake file and sources from Dummy:
sed -e "s/Dummy/$MODNAME/g" $MODDIR/Dummy/CMakeLists.txt > $MODDIR/$MODNAME/CMakeLists.txt

# Copy over the README, setting current git username/email as author
# If this fails, use system username and hostname
MYNAME=$(git config user.name)
MYMAIL=$(git config user.email)
if [ -z "$MYNAME" ]; then
    MYNAME=$(whoami)
fi
if [ -z "$MYMAIL" ]; then
    MYMAIL=$(hostname)
fi

if [ "$type" == 1 ]; then
    MODTYPE="\*GLOBAL\*"
fi
if [ "$type" == 2 ]; then
    MODTYPE="*DETECTOR* **Detector Type**: *<add types here>*"
fi
if [ "$type" == 3 ]; then
    MODTYPE="\*DUT\*"
fi

sed -e "s/Dummy/$MODNAME/g" \
    -e "s/\*NAME\*/$MYNAME/g" \
    -e "s/\*EMAIL\*/$MYMAIL/g" \
    -e "s/\*GLOBAL\*/$MODTYPE/g" \
    -e "s/Functional/Immature/g" \
    "$MODDIR/Dummy/README.md" > "$MODDIR/$MODNAME/README.md"

# Copy over source code skeleton:
sed -e "s/Dummy/$MODNAME/g" "$MODDIR/Dummy/Dummy.h" > "$MODDIR/$MODNAME/${MODNAME}.h"
sed -e "s/Dummy/$MODNAME/g" "$MODDIR/Dummy/Dummy.cpp" > "$MODDIR/$MODNAME/${MODNAME}.cpp"

# Options for sed vary slightly between mac and linux
opt=-i
platform=`uname`
if [ "$platform" == "Darwin" ]; then opt="-i \"\""; fi

# Change to detector or dut module type if necessary:
if [ "$type" != 1 ]; then

    # Prepare sed commands to change to per detector module
    # Change module type in CMakeLists
    if [ "$type" == 3 ]; then
        command="sed $opt 's/_GLOBAL_/_DUT_/g' $MODDIR/$MODNAME/CMakeLists.txt"
    else
        command="sed $opt 's/_GLOBAL_/_DETECTOR_/g' $MODDIR/$MODNAME/CMakeLists.txt"
    fi
    eval $command

    # Add aux exclusion for detector module if selected:
    if [ "$exclude_aux" == 1  ]; then
        command="sed $opt '/_DETECTOR_/ a CORRYVRECKAN_EXCLUDE_AUX(\${MODULE_NAME})' $MODDIR/$MODNAME/CMakeLists.txt"
        eval $command
    fi

    # Add inclusion of passive detectors for detector module if selected:
    if [ "$exclude_pass" == 0  ]; then
        command="sed $opt '/_DETECTOR_/ a CORRYVRECKAN_INCLUDE_PASSIVE(\${MODULE_NAME})' $MODDIR/$MODNAME/CMakeLists.txt"
        eval $command
    fi

    # Change header file
    command="sed ${opt} \
    -e 's/param detectors.*/param detector Pointer to the detector for this module instance/g' \
    -e 's/std::vector<std::shared_ptr<Detector>> detectors/std::shared_ptr<Detector> detector/g' \
    $MODDIR/$MODNAME/${MODNAME}.h"
    eval $command
    # Change implementation file
    command="sed ${opt} \
    -e 's/std::vector<std::shared_ptr<Detector>> detectors/std::shared_ptr<Detector> detector/g' \
    -e 's/move(detectors)/move\(detector\)/g' \
    $MODDIR/$MODNAME/${MODNAME}.cpp"
    eval $command
fi

# Print a summary of the module created:
FINALPATH=`realpath $MODDIR/$MODNAME`
echo "Name:   $MODNAME"
echo "Author: $MYNAME ($MYMAIL)"
echo "Path:   $FINALPATH"
echo
echo "Re-run CMake in order to build your new module."
