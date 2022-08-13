#!/bin/bash

mkdir ../pandelos_output/example1

read -p "Scrivi i numeri " -a input

IFS=', ' read -a ranges <<< "${input[@]}"
for range in "${ranges[@]}"; do
        IFS=- read start end <<< "$range"
        [ -z "$start" ] && continue
        [ -z "$end" ] && end=$start
        for (( i=start ; i <= end ; i++ )); do
                bash pandelos.sh 'example1/mycoplasma_'$i'.faa' 'example1/output_mycoplasma_'$i'.net' 1 'example1/log_mycoplasma_'$i'.txt' 'example1/mycoplasma_'$i'.clus'
        done
done