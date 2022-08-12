#!/bin/bash

g++ -w -std=c++17 -fopenmp -O3 src/cpp/main.cpp -o bin/pandelos_plus_plus.out

echo "Finished build"

while true; do
    read -p "Do you want to download the sample files from the internet? " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

[ ! -d "../pandelos_data" ] && mkdir ../pandelos_data

[ ! -d "../pandelos_output" ] && mkdir ../pandelos_output

mkdir ../pandelos_data/mycoplasma_panprova

wget https://pangenes.s3.eu-south-1.amazonaws.com/mycoplasma_panprova/mycoplasma_panprova_16.zip -O ../pandelos_data/mycoplasma_panprova/mycoplasma_panprova_16.zip

#wget http://localhost:8888/mycoplasma_panprova_16.zip

unzip ../pandelos_data/mycoplasma_panprova/mycoplasma_panprova_16.zip -d ../pandelos_data/mycoplasma_panprova/ && rm ../pandelos_data/mycoplasma_panprova/mycoplasma_panprova_16.zip

mkdir ../pandelos_data/escherichia_panprova

wget -O ../pandelos_data/escherichia_panprova/escherichia_panprova_16.zip https://pangenes.s3.eu-south-1.amazonaws.com/escherichia_panprova/escherichia_panprova_16.zip

unzip ../pandelos_data/escherichia_panprova/escherichia_panprova_16.zip -d ../pandelos_data/escherichia_panprova/ && rm ../pandelos_data/escherichia_panprova/escherichia_panprova_16.zip

mkdir ../pandelos_data/example1 && mkdir ../pandelos_data/example2 && mkdir ../pandelos_data/example3 && mkdir ../pandelos_data/example4

cd ..

current_path=$(pwd)

#
# It deletes unnecessary files and transforms PANPROVA's output into .gbk files that you will need later
# 
echo "PANPROVA to gbk"

bash PanDelos-plusplus/script/panprova2gbk.sh $current_path'/pandelos_data/mycoplasma_panprova'

bash PanDelos-plusplus/script/panprova2gbk.sh $current_path'/pandelos_data/escherichia_panprova'

#
# Transform .gbk files into .faa files composed of amino acid sequences
# $1 = absolute path where the .gbk files are present
# $2 = absolute path and prefix of the .faa file to generate
#
echo "Transform .gbk files into .faa files"

bash PanDelos-plusplus/script/gbk2faa.sh $current_path'/pandelos_data/mycoplasma_panprova/' $current_path'/pandelos_data/example3/mycoplasma_panprova_'

bash PanDelos-plusplus/script/gbk2faa.sh $current_path'/pandelos_data/escherichia_panprova/' $current_path'/pandelos_data/example4/escherichia_panprova_'

#
# Transform the output of PANPROVA into a faa file consisting of nucleotide sequences
# $1 = input of the files generated by PANPROVA
# $2 = common prefix of all genomes
# $3 = path and common prefix of all .faa files to be generated
#
while true; do
    read -p "Do you want to generate faa files with corresponding nucleotides? It may take several minutes" yn
    case $yn in
        [Yy]* ) bash PanDelos-plusplus/script/panprova2nucleotides.sh $current_path'/pandelos_data/mycoplasma_panprova/' 'genome_' $current_path'/pandelos_data/example1/mycoplasma_';
                bash PanDelos-plusplus/script/panprova2nucleotides.sh $current_path'/pandelos_data/escherichia_panprova/' 'genome_' $current_path'/pandelos_data/example2/escherichia_';
                break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done