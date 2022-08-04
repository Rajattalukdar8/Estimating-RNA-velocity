# !/bin/bash


GTF=home/ansari20336/anaconda3/bin/hg19.ncbiRefSeq.gtf

FILES=/storage/vibhor/Mtech/Ariba/BamFiles/*.bam
OUTPUT=/storage/vibhor/Mtech/Ariba/Loom


for f in $FILES
do
   echo = $f
  ./velocyto run -c -U -o $OUTPUT -m hg19_rmsk.gtf $f /home/ansari20336/anaconda3/bin/hg19.ncbiRefSeq.gtf
done


echo "done!"



