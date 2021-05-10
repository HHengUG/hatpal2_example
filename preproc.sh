#!/bin/bash
echo $1 
# input BAM
echo $2
# output directory
echo $3
# project name

filtered=$2$3"_F.bam"
sorted=$2$3"_S.bam"
psam=$2$3"_P3S.sam"
nsam=$2$3"_N3S.sam"
header=$2$3"_HEADER"
sbam=$2$3"_3S.bam"
ssorted=$2$3"_S3S.bam"

samtools view -@ 2 -F 256 -b -h $1 > $filtered 
samtools sort $filtered > $sorted  
samtools index -@ 2 $sorted  
samtools view -F 16 $sorted | awk '{if ($6~/([0-9][0-9]|[3-9])S$/ ) print }' > $psam 
samtools view -f 16 $sorted | awk '{if ($6~/^([0-9][0-9]|[3-9])S/ ) print }' > $nsam 
samtools view -H $sorted > $header 
cat $header $psam $nsam | samtools view -b > $sbam
samtools sort $sbam > $ssorted 
samtools index -@ 2 $ssorted 

echo $ssorted
