# to use STAR to align the reads to the reference genome
if [ $# -eq 0 ]
then
	echo "sh batch_fastp_process.sh target_file post_fix"
    echo "post_fix can be fastq  fq  fq.gz fastq.gz, Default: fastq"
	exit 0
fi



if [ $# -eq 2 ]
then
	post_fix=$2
else
	post_fix="fastq"
fi

if [ ! -d ./fastp ]
then
    mkdir ./fastp
fi

for file_prefix in `cat $1`
do
    fastp \
    -i ${file_prefix}_1.${post_fix} -I ${file_prefix}_2.${post_fix} \
    -o ./fastp/${file_prefix}_clean_1.fastq -O ./fastp/${file_prefix}_clean_2.fastq \
    --detect_adapter_for_pe \
    -q 20 \
    -u 70 \
    -n 5 \
    -l 36 \
    -h report-${file_prefix}.html \
    -j report-${file_prefix}.json
done
