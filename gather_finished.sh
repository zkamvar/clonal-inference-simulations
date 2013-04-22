#!/bin/sh

mkdir finals;
for i in clone*/;
	do 
		for j in $i/clone*gen_10.0k*dat;
			do
				cp $j finals;
		done
done
        
