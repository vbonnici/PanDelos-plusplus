#!/bin/bash

absolute_path_in=$1
output_file=$2

for i in $(seq 2 2 16)
do
    python3 PanDelos-plusplus/script/gbk2faa.py $absolute_path_in $output_file$i'.faa' $i
done