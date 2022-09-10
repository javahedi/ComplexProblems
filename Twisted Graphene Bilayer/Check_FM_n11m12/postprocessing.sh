#!/bin/bash


dir=$(pwd)
cd $dir/job
rm *.dat

for i in {0..120}    
    	do
                cd U$i
                #rm HK
		fileIn="Magnetization"
		while IFS= read -r line
			do
        			# display $line or do somthing with $line
				#printf '%s\n' "$line"
               			echo "$line"  >> ../Mz_75.dat
				#echo $pwd
                done <"$fileIn"

               
		fileIn="E_Gap_u.dat"
		while IFS= read -r line
			do
        			# display $line or do somthing with $line
				#printf '%s\n' "$line"
               			echo "$line"  >> ../Gap_75.dat
				#echo $pwd
                done <"$fileIn"                

                
		fileIn="E_GState.dat"
		while IFS= read -r line
			do
        			# display $line or do somthing with $line
				#printf '%s\n' "$line"
               			echo "$line"  >> ../GroundEn_75.dat
				#echo $pwd
                done <"$fileIn"


                cd ..

done
