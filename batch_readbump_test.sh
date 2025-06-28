if [ $# -lt 4 ]
then
	echo "sh batch_readbump_test.sh target_dir target_file observing_window_file lib_type"
    # readsends -w observing_window.tsv -t sam -f condition1.sam -c test
    echo "target_dir: is the folder where the bam files are located, ending with /"
    echo "target_file: is the file the contains the prefix of the bam files"
    echo "observing_window_file: is the file that contains the range to test the read ends"
    echo "lib_type: should be 1 or 2 or 0"
	exit 0
fi

target_dir=$1
lib_type=$4
observing_window=$3
target_file=$2
if ! [ -f $observing_window ]
then
    echo " ${observing_window} not found"
    exit 0
fi 

for name in `cat ${target_file}`
do
    # sorted_${name}.bam
    sam_file="${target_dir}${name}_subset.sam"    
    if [ -f ${sam_file} ]
    then        
        echo "+testing ${sam_file} with ${observing_window} to bumps_details_${name}.tsv "
        # readsends -w ${observing_window} -t sam -f ${sam_file} -c test -r -l ${lib_type} >bumps_${name}.tsv
        readsends -w ${observing_window} -t sam -j 60 -u 70 -f ${sam_file} -c test  -l ${lib_type} >bumps_details_${name}.tsv 
    else
         echo "?? failed ${sam_file} with ${observing_window} to bumps_details_${name}.tsv"
    fi   
done
