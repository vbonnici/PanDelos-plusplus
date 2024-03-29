#!/bin/bash


echo "================================================================================"
echo "PanDelos plus plus example4"
echo -e
echo "In this example you can analyze up to 16 escherichia coli genomes with amino acid sequences"
echo "You will need to enter one or more numbers of genomes you intend to analyze. Examples:"
echo -e
echo "I enter the numbers: 4, 6"
echo "Run 2 separate tests with escherichia_4.faa (genome 0 to genome 3) and escherichia_6.faa (analogue)"
echo -e
echo "I enter the numbers: 2-4, 6, 8"
echo "Run 5 separate tests taking into account the files ranging from escherichia_2.faa to escherichia_4.faa and the 2 separate files escherichia_6.faa and escherichia_8.faa"
echo "================================================================================"

[ ! -d "../pandelos_output/example4" ] && mkdir ../pandelos_output/example4

read -p "Write the numbers of the genomes you want to analyze " -a input

IFS=', ' read -a ranges <<< "${input[@]}"
for range in "${ranges[@]}"; do
        IFS=- read start end <<< "$range"
        [ -z "$start" ] && continue
        [ -z "$end" ] && end=$start
        for (( i=start ; i <= end ; i++ )); do
             [ $((i%2)) -eq 0 ] && bash pandelos.sh 'example4/escherichia_'$i'.faa' 'example4/output_escherichia_'$i'.net' 0 'example4/log_escherichia_'$i'.txt' 'example4/escherichia_'$i'.clus' 'example4/escherichia_coco_'$i'.txt'
        done
done