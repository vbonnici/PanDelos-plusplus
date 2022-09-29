#!/bin/bash


echo "================================================================================"
echo "PanDelos plus plus example1"
echo -e
echo "In this example you can analyze up to 16 mycoplasma genitalium genomes with nucleotide sequences"
echo "You will need to enter one or more numbers of genomes you intend to analyze. Examples:"
echo -e
echo "I enter the numbers: 4, 6"
echo "Run 2 separate tests with mycoplasma_4.faa (genome 0 to genome 3) and mycoplasma_6.faa (analogue)"
echo -e
echo "I enter the numbers: 2-4, 6, 8"
echo "Run 5 separate tests taking into account the files ranging from mycoplasma_2.faa to mycoplasma_4.faa and the 2 separate files mycoplasma_6.faa and mycoplasma_8.faa"
echo "================================================================================"

[ ! -d "../pandelos_output/example1" ] && mkdir ../pandelos_output/example1

read -p "Write the numbers of the genomes you want to analyze " -a input

IFS=', ' read -a ranges <<< "${input[@]}"
for range in "${ranges[@]}"; do
        IFS=- read start end <<< "$range"
        [ -z "$start" ] && continue
        [ -z "$end" ] && end=$start
        for (( i=start ; i <= end ; i++ )); do
             [ $((i%2)) -eq 0 ] && bash pandelos.sh 'example1/mycoplasma_'$i'.faa' 'example1/output_mycoplasma_'$i'.net' 1 'example1/log_mycoplasma_'$i'.txt' 'example1/mycoplasma_'$i'.clus' 'example1/mycoplasma_coco_'$i'.txt'
        done
done