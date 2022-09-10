#!/bin/bash

dir=$(pwd)
cd $dir

for d in */ ; do
    #echo "$d"
    cd $d
    for m in */ ; do
        #echo "$m"
        cd $m
            rm -r __pycache__ lanczos sub* *.py N*
        cd ..
    done
    cd ..
done
