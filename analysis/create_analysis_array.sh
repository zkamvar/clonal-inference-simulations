#!/bin/env bash
#
# This script will create a text file containing commands to run an array
# using the unzip_and_analyze script.
if [ $# -lt 3 ]; then
	echo 
	echo "Create a file of commands for analysis with SGE_Array"
	echo
	echo "Usage:"
	echo
	echo "    bash create_analysis_array.sh <path> <script> <outfile>"
	echo
	echo "    <path>    - path to zipped feather files"
	echo "    <script>  - script with arguments to run on each file once unzipped"
	echo "    <outfile> - text file to place all the commands into"
	echo
	exit
fi


zipDir=$1
theScript=$2
theFile=$3
controlScript="bash unzip_and_analyze.sh "

if [ -e "$theFile" ]; then
	if [ -s "$theFile" ]; then
		echo "File "$theFile" exists and is not empty. Exiting."
		exit
	fi
else 
	touch $theFile
fi

for i in $(ls -1 -d "$zipDir"/*zip)
do
	echo $controlScript$i" \""$theScript"\"" >> $theFile
done

nlines=$(cat "$theFile" | wc -l)

echo "successfully gathered "$nlines" line(s)."