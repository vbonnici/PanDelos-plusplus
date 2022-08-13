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
echo -e
echo -e

input_file=$1
output_file=$2
sequences_type=$3
log_file=$4
clus_file=$5

echo "Clustering: "
date
bin/pandelos_plus_plus.out -f '../pandelos_data/'$input_file -o '../pandelos_output/'$output_file -s $sequences_type -l '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
date

echo -e
echo -e

echo "De-clustering: "
date
echo "Start netclu_ng_parallel " >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
python3 src/python/netclu_ng_parallel.py '../pandelos_data/'$input_file '../pandelos_output/'$output_file >> '../pandelos_output/'$log_file
echo "End netclu_ng_parallel " >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
date

echo -e
echo -e

echo "Writing gene gene families in file clus: "
date
grep "F{ " '../pandelos_output/'$log_file | sed s/F{\ //g | sed s/}//g | sed s/\ \;//g | sort | uniq >'../pandelos_output/'$clus_file
date

echo -e
echo -e

echo "Quality: "
date
echo "Quality: " >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
python3 src/python/quality.py '../pandelos_data/'$input_file '../pandelos_output/'$clus_file >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
echo -e >> '../pandelos_output/'$log_file
date

echo -e
echo -e

echo "Finish!"
echo -e
echo -e

if (( $SECONDS > 3600 )) ; then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo -e >> '../pandelos_output/'$log_file
    echo -e >> '../pandelos_output/'$log_file
    echo "TOTAL TIME: $hours hour(s), $minutes minute(s) and $seconds second(s)" >> '../pandelos_output/'$log_file
elif (( $SECONDS > 60 )) ; then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo -e >> '../pandelos_output/'$log_file
    echo -e >> '../pandelos_output/'$log_file
    echo "TOTAL TIME: $minutes minute(s) and $seconds second(s)" >> '../pandelos_output/'$log_file
else
    echo -e >> '../pandelos_output/'$log_file
    echo -e >> '../pandelos_output/'$log_file
    echo "TOTAL TIME: $SECONDS seconds" >> '../pandelos_output/'$log_file
fi

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