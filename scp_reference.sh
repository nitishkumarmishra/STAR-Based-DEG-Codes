#!/bin/bash

RPROF=/research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish
SCRATCH=/research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish/scratch
TMPPROJDIR=$SCRATCH/tmpprojdir


TMPPROJDIRrsem="$TMPPROJDIR/star-ref"
mkdir  -p  $TMPPROJDIRrsem
TMPPROJDIRstar="$TMPPROJDIR/star-ref"
mkdir  -p  $TMPPROJDIRstar


TMPWF="$TMPPROJDIR/tmp.work-albqsr.$BASHPID"
TMPPROJDIRdata="$TMPWF.tmp-data-tmp"
mkdir  -p  $TMPPROJDIRdata



GENOMEDIR="$RPROF/data/mouse/ensembl"
GENOMEFA="Mus_musculus.GRCm39.dna.primary_assembly.rdna.fa.gz"
GENOMEFAtag=$(basename $GENOMEFA  .gz)
GTFtag="Mus_musculus.GRCm39.104.rdna_rn18s.gtf"


mkdir  -p  $TMPPROJDIRdata
rsync  -tvhL  $GENOMEDIR/$GENOMEFAtag.gz  $TMPPROJDIRdata
gunzip  $TMPPROJDIRdata/$GENOMEFAtag.gz

GTFbase=$(basename  $GTFtag)
rsync  -vhL  $GENOMEDIR/$GTFtag.gz  $TMPPROJDIRdata
gunzip  -f $TMPPROJDIRdata/$GTFtag.gz

STARDIR="$RPROF/data/STAR"
STARTAG="STAR.index.$GENOMEFAtag--$GTFtag"

RSEMDIR="$RPROF/data/rsem-transcriptome"
RSEMTAG="$GENOMEFAtag--$GTFtag.rsem.transcriptome"



############ RSEM index
echo -e "rsem index...\n"
if [ -d $TMPPROJDIRrsem ] ; then
rsync  -tvh  $RSEMDIR/$RSEMTAG.tar.gz     $TMPPROJDIRrsem

echo -e "\tunpack index...\n"

tar  zxv  -C $TMPPROJDIRrsem  -f $TMPPROJDIRrsem/$RSEMTAG.tar.gz
fi



############ STAR index
echo -e "STAR index...\n"
if [ -d $TMPPROJDIRstar ]; then
rsync  -tvh  $STARDIR/$STARTAG.tar.gz   $TMPPROJDIRstar

echo -e "\tunpack index...\n"

tar  zxv  -C $TMPPROJDIRstar  -f $TMPPROJDIRstar/$STARTAG.tar.gz
fi

