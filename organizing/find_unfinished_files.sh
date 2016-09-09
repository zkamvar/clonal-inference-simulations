#!/bin/env bash
#
# 2016-09-09
# This script will travel into a directory, look through the error files for a
# specific error pattern, and return a new file to use for re-running the 
# analysis.
# Use:
# bash zhian_simulations/organizing/find_unfinished_files.sh <path> <error message> <outfile>

dataDir=$1
error=$2
theFile=$3

if [ -e "$theFile" ]; then
	if [ -s "$theFile" ]; then
		echo "File "$theFile" exists and is not empty. Exiting."
		exit
	fi
else 
	touch $theFile
fi

cmd=$dataDir"commands.txt"

# The prefix is obtained by cutting the trailing slash and cutting any preceding
# directory structure.
prefix=$(echo "$dataDir" | rev | cut -c 2- | cut -d "/" -f -1 | rev)
efiles=$dataDir$prefix".e*"

ERRORS=$(grep -il $error $efiles)

for e in $ERRORS
do
	echo $e
	errno=$(echo $e | rev | cut -d "." -f -1 | rev) # getting the line number
	sed -n $errno"p" < $cmd >> $theFile             # extract line from $cmd
done
