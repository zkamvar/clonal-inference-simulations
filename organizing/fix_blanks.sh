#!/usr/bin env sh

FOLDERS=pt*/

for d in $FOLDERS
do
	echo $d
	cd $d
	for f in *pop
	do
		new=$(echo $f | tr ' ' '0')
		mv "-T" "$f" "$new"
	done
	cd ..
done
