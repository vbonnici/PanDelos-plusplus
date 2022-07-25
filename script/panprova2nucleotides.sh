#!/bin/bash

absolute_path_in=$1
prefix=$2
output_file=$3

for i in $(seq 2 2 16)
do
    a=$(($i -1))
    for j in $(seq 0 $a)
	do
	  	python3 pandelos-plus-plus/script/panprova2nucleotides.py $absolute_path_in$prefix$j'.fna' $absolute_path_in$prefix$j'.gff' $output_file$i'.faa'
	done
done