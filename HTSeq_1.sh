#!/bin/sh
start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

#sample=samples.txt
#REF=/work/unmc_gudalab/nitish123/STAR_hg38/STAR_hg38_V33/
#GTF=/work/unmc_gudalab/nitish123/STAR_hg38/gencode.v33.annotation.gtf
outpath='Dataset1_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}
for sample in $(ls Dataset1 |grep "fastq"| sed s/[12].fastq.gz// | sort -u);
do echo
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
STAR \
--readFilesIn Dataset1/${sample}1.fastq.gz Dataset1/${sample}2.fastq.gz \
--runThreadN 15 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignIntronMin 20 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimOutType WithinBAM \
--chimSegmentMin 15 \
--genomeDir /work/unmc_gudalab/nitish123/STAR_hg38/STAR_hg38_V33/ \
--genomeLoad NoSharedMemory \
--limitSjdbInsertNsj 1200000 \
--outFilterIntronMotifs None \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterType BySJout \
--outSAMattributes NH HI AS nM NM ch \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat \
--outFileNamePrefix ${outpath}/${sample}/${sample}_ \
--twopassMode Basic;

done
end=`date +%s`
runtime=$((end-start))
echo $runtime

