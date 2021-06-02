rm -rf stat_mapping.txt
rm -rf ~/temp/
export TMPDIR=~/temp/
while read sample;
do
mapped=$(samtools view -F 0x4 alignment/${sample}.sort.bam | cut -f 1 | sort | uniq | wc -l)
total=$(samtools view alignment/${sample}.bam | wc -l)
echo -e $sample"\t"$((total/2))"\t"$mapped >> stat_mapping.txt
done < data_ID.txt

