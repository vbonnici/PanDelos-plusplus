echo "################################################################################"
echo "# PanDelos plus plus:                                                          #"
echo "# a dictionary-based method for pan-genome content discovery.                  #"
echo "#------------------------------------------------------------------------------#"
echo "# Bonnici et. al                                                               #"
echo "# Giandonato Inverso                                                           #"
echo "#==============================================================================#"
echo "# This software is under MIT license!                                          #"
echo "# Please visit https://github.com/vbonnici/PanDelos-plusplus                   #"
echo "################################################################################"
echo "\n"
echo "\n"

input_file=$1
output_file=$2
sequences_type=$3
log_file=$4
clus_file=$5

echo "clustering ..."
date
bin/pandelos_plus_plus.out -f '../pandelos_data/'$input_file -o '../pandelos_output/'$output_file -s $sequences_type -l '../pandelos_output/'$log_file
date
echo "de-clustering ..."
python3 src/python/netclu_ng_parallel.py '../pandelos_data/'$input_file '../pandelos_output/'$output_file >> '../pandelos_output/'$log_file

echo "writing gene gene families in file clus ..."
date
grep "F{ " '../pandelos_output/'$log_file | sed s/F{\ //g | sed s/}//g | sed s/\ \;//g | sort | uniq >'../pandelos_output/'$clus_file

python3 src/python/quality.py '../../../pandelos_data/'$input_file '../../../pandelos_output/'$clus_file

echo "Finish!"
