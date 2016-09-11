#!/bin/env bash
#
# This script is used to unzip a file to the temporary directory and then
# execute a given script
#
# Update 2016-09-10:
# 
# This script will clean up after itself now. 
# The {} || {} style of a try catch came from this answer on stackoverflow:
# http://stackoverflow.com/a/22010339/2752888
#
if [ $# -lt 2 ]; then
	echo
	echo "Unzip a file in the temporary directory of a node and analyze the contents with a script"
	echo
	echo "Usage:"
	echo
	echo "    bash unzip_and_analyze.sh <file> <script>"
	echo
	echo "    <file>   - zipped feather format file to analyze"
	echo "    <script> - script with arguments to run on each file once unzipped"
	echo
	exit
fi

theFile=$1
theScript=$2

# If there are, create a new directory and assign that to TMPDIR
if [ "$(ls -A $TMPDIR)" ]; then
	mkdir $TMPDIR/znk
	TMPDIR=$TMPDIR/znk
fi

{
	# Unzip the file into TMPDIR
	unzip -j -d $TMPDIR $theFile
} || {
	echo "removing "$TMPDIR
	rm -rf $TMPDIR
	exit
}

# Create the command by concatenating the results of the unzipping
# and the script to run. Note, for multiple arguments within a
# script, quoting is encouraged (e.g. "ls -l")
cmd=$theScript" "$(ls -d $TMPDIR/*)
echo $cmd
{
	eval $cmd
} || {
	echo "removing "$TMPDIR
	rm -rf $TMPDIR
	exit
}

echo "removing "$TMPDIR
rm -rf $TMPDIR
exit