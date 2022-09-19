#!/bin/bash
#TQ="WT_1/ERR458493.fastq.gz,WT_1/ERR458494.fastq.gz,WT_1/ERR458495.fastq.gz,WT_1/ERR458496.fastq.gz"
TQ=$(find -name "*L00*R1*.fastq.gz" | sort | tr "\n" ","|sed "s/\.\///g"|sed 's/,$//')
echo $TQ
