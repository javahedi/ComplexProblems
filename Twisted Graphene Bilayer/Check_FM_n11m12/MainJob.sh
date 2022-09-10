#!/bin/bash

dir=$(pwd)
a=0.2  #step of U potential
s=0.5

cd $dir
mkdir job
cd job

for i in {1..4}
    do       
       
       mkdir U$i
       cd U$i          
          echo $s >> Hubard_U.dat
          cp -a ../../program/. .        
          qsub subrout.job
            
       cd ..
      
    s="$( bc <<<"$s + $a" )"       
done
