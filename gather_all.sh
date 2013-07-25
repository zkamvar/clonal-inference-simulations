#!/bin/sh

mkdir finals;
for i in clone*/;
	do 
		for j in $i/clone*dat;
			do
				cp $j finals;
		done
done
        
