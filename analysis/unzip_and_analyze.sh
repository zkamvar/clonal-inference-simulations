#!/bin/env bash
#
# This script is used to unzip a file to the temporary directory and then
# execute a given script

theFile=$1
theScript=$2

# If there are, create a new directory and assign that to TMPDIR
if [ "$(ls -A $TMPDIR)" ]; then
	mkdir $TMPDIR/znk
	TMPDIR=$TMPDIR/znk
fi

# Unzip the file into TMPDIR
unzip -j -d $TMPDIR $theFile

# Create the command by concatenating the results of the unzipping
# and the script to run. Note, for multiple arguments within a
# script, quoting is encouraged (e.g. "ls -l")
cmd=$theScript" "$(ls -d $TMPDIR/*)
echo $cmd
eval $cmd
