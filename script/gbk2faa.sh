#!/bin/bash

absolute_path_in=$1
output_file=$2

for i in $(seq 2 2 16)
do
    a=$(($i -1))
    for j in $(seq 0 $a)
    do
        python3 pandelos-plus-plus/script/gbk2faa.py $absolute_path_in $output_file$i'.faa'
    done
done