#!/bin/bash

dir=$1
cd $dir

rm *.fna;
rm *.gff;

find . -type f -name "*.gbff" -exec rename 's/\.gbff$/.gbk/' '{}' \;

for filename in ./*.gbk; do
	result="${filename}.gbk"
	#sed "2 i ACCESSION ${filename:2}" $filename > $result
done

for filename in *
do
    x=`echo $filename | sed 's/\..*\./\./g'`
    mv $filename $x
done
