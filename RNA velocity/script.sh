# !/bin/bash


index=/storage/vibhor/Mtech/Ariba/genome_dir/genomehg19_index

FILES=/storage/vibhor/Mtech/Ariba/FastqFiles/*_1.fastq
OUTPUT=/storage/vibhor/Mtech/Ariba/BamFiles

for f in $FILES
do

 #  echo = $f
   base=$(basename $f.)

#    echo = $base
    echo = $f ${f%_1.fastq}_2.fastq

   ./STAR --runThreadN 10 --genomeDir $index --readFilesIn  $f ${f%_1.fastq}_2.fastq\
          --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --outFileNamePrefix $OUTPUT/$base
done


echo "done!"
