# to use STAR to align the reads to the reference genome
if [ $# -lt 2 ]
then
	echo "sh batch_align_sort_human.sh prefix_name_file postfix read_count"
	echo "postfix should be like fastq fastq.gz fq fq.gz"
	echo "read_count: 1 or 2"
	exit 0
fi
# change the genome dir below
genome_dir="/home/jhy/NGS_data/GRch38/star_db"

postfix=$2
data_reader=""
if [[ postfix =~ /gz$/ ]]
then
	data_reader="zcat"
fi
	

if [ $# -eq 3 ]
then
	read_count=$3
else
	read_count=2
	# default is paired reads
fi



for name in `cat $1`
do
	if [ ${read_count} -eq 1 ]
	then
		file_name_prefix1="${name}.${postfix}"
		file_name_prefix2=""
	else
		file_name_prefix1="${name}_1.${postfix}"
		file_name_prefix2="${name}_2.${postfix}"
	fi

	if [ -f "${file_name_prefix1}" ]; then
		echo "Aligning ${file_name_prefix1} ${file_name_prefix2} using STAR:"
		STAR \
		--runThreadN 60 \
		--genomeDir ${genome_dir} \
		--readFilesIn "${file_name_prefix1}" "${file_name_prefix2}" \
		--outFileNamePrefix "result_${name}_" \
		--limitBAMsortRAM 1000000000000 \
		--genomeLoad LoadAndKeep \
		--outSAMtype BAM Unsorted \
		# --readFilesCommand ${data_reader}

		echo "Sorting sorted_${name}.bam using samtools"
		samtools sort -@70 -o "sorted_${name}.bam" -T "/tmp/samtool_${name}" -O BAM "result_${name}_Aligned.out.bam"
		
		echo "Deleting result_${name}_Aligned.out.bam"
		rm  result_${name}_Aligned.out.bam

		echo "Indexing bam file:sorted_${name}.bam"
		samtools index -b sorted_${name}.bam
	else
		echo "${file_name_prefix1} does not exists"
	fi
done
