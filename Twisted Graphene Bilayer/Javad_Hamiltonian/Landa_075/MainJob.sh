#!/bin/bash

dir=$(pwd)
a=0.04  #step of U potential
s=0.0

cd $dir
mkdir job_nk400
cd job_nk400

for i in {0..25}
    do        
       #rm -r U$i
       mkdir U$i
       cd U$i          
          echo $s >> Hubard_U.dat
       cp -a ../../program/. . 
       #qsub subrout.uge
       qsub subrout.job
       #rm core.*       
       cd ..
       
       
    #s=`echo $s + $a | bc`
    #s="$s + $a" | bc
    s="$( bc <<<"$s + $a" )"       
done
