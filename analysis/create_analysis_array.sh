#!/bin/env bash
#
# This script will create a text file containing commands to run an array
# using the unzip_and_analyze script.

zipDir=$1
theScript=$2
theFile=$3
controlScript="bash zhian_simulations/analysis/unzip_and_analyze.sh "

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