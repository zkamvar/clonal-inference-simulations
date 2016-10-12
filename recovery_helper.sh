#!/bin/bash

if [ $# -eq 0 ]; then
  echo
  echo "Rerun unused simulations"
  echo
  echo "Usage:"
  echo
  echo "    bash recovery_helper.sh <failed simulation text file>"
  echo
  exit
fi

# http://stackoverflow.com/a/10929511
foo=( )

counter=0
while IFS='' read -r line || [[ -n "$line" ]]; do
    # echo "Text read from file: $line; counter $counter"
    foo[$counter]=$line
    let counter=counter+1
done < "$1" 

PIDS=( )
counter=0
# http://stackoverflow.com/a/6723516
for i in "${!foo[@]}"; do 
  if ! ((i % 8)); then
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
