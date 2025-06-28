if [ $# -lt 3 ]
then
	echo "sh batch_extract_reads.sh target_dir target_file input_bed"
    echo "target_dir: is the folder where the bam files are located, ending with /"
    echo "target_file: is the file the contains the prefix of the bam files"
    echo "input_bed: is the bed file that contains the extracting range of the bam files"
	exit 0
fi
target_dir=$1
for name in `cat $2`
do
    bam_file="${target_dir}${name}.bam"
    input_bed=$3
    if [ -f ${bam_file} ]
    then
        if [ -f $input_bed ] 
        then
            echo "+extracting ${bam_file} with ${input_bed} to ${name}_subset.bam"
            samtools view -hb -L ${input_bed} ${bam_file} > ${name}_subset.bam
            echo " converting ${name}_subset.bam to ${name}_subset.sam"
            samtools view -h -o ${name}_subset.sam  ${name}_subset.bam
            # echo " converting ${name}_subset.bam to ${name}_subset.bed"
            # bedtools bamtobed -i  ${name}_subset.bam > ${name}_subset.bed
            echo " $bam_file done" 
        else
            echo " ${input_bed} not found"
        fi
    else
        echo " ${bam_file} not found"
    fi 
done
