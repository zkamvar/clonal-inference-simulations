#!/bin/bash

if [ $# -lt 2 ]; then
	echo
	echo "Rerun unused simulations"
	echo
	echo "Usage:"
	echo
	echo "		bash recovery_helper.sh <number of concurrent processes> <failed simulation text file>"
	echo
	exit
fi

# Reading in file line by line
# http://stackoverflow.com/a/10929511
foo=( )

counter=0
while IFS='' read -r line || [[ -n "$line" ]]; do
	# echo "Text read from file: $line; counter $counter"
	foo[$counter]=$line
	let counter=counter+1
done < "$2" 

PIDS=( )
counter=0
# Looping over an array with index
# http://stackoverflow.com/a/6723516
for i in "${!foo[@]}"; do 
	if ! ((i % $1)); then
		# For the PID waiting
		# http://stackoverflow.com/a/356154
		let counter=0
		for pid in "${!PIDS[@]}"; do
			printf "waiting on ${PID[$pid]}\n"
			wait ${PID[$pid]}
		done
	fi
	printf "%s\t%s\n" "$i" "${foo[$i]}"
	eval ${foo[$i]} & 
	PIDS[$counter]=$!
	let counter=counter+1
done
wait
