#!/bin/bash

source bashrc_rdna_variant_calling

echo -e "\n\n---------------- generate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh -----------------\n"
echo -e "This script will generate a Genome Index for mRNA for STAR\n"
echo -e "\t\tcommand:: sh generate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh GENCODEmm39 -p 10\n"


while  test  $# -gt 0;  do
        case "$1"  in
                -h|--help|-help)
                        echo -e "\n\ngenerate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh\n"
                        echo -e "\n\n"
                        exit  0
			;;
                -p)
                        # threads
                        shift
                        if test $# -gt 0; then
                                NUMTHREADS="$1"
                        else
                                echo -e "\n\nERROR! number of threads not specified.\n\n"
                                exit 1
                        fi
                        shift
                        ;;
  		*) break;;
        esac
done

GENOMETAG="$1"

############################ check parameters
echo -e "\t\t\tGENOMETAG=[$GENOMETAG]\n"
############################ check parameters
echo -e "\n\n\n"
echo -e "GENOMETAG=[$GENOMETAG]"
echo -e "NUMTHREADS=[$NUMTHREADS]"
echo -e "\n\n\n"


SPECIES=$GENOMETAG


case $SPECIES in
  hg38|h38)
    GENOMEDIR="$RPROF/data/human/ensembl"
    GENOMEFA="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GENOMEFAtag=$(basename $GENOMEFA  .gz)
    GTFtag="Homo_sapiens.GRCh38.97.rdna_rn18s.gtf"
    ;;
  GENCODEhg38|GENCODEh38)
    GENOMEDIR="$RPROF/data/human/gencode"
    GENOMEFA="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GENOMEFAtag=$(basename $GENOMEFA  .gz)
    GTFtag="Homo_sapiens.GRCh38.97.rdna_rn18s.gtf"
    ;;
  mm10|m38)
    GENOMEDIR="$RPROF/data/mouse/mm10"
    GENOMEFA="Mus_musculus.GRCm38.dna.primary_assembly.rdna.fa.gz"
    GENOMEFAtag=$(basename $GENOMEFA  .gz)
    GTFtag="Mus_musculus.GRCm38.97.rdna_rn18s.gtf"
    ;;
  mm39|m39)
    GENOMEDIR="$RPROF/data/mouse/ensembl"
    GENOMEFA="Mus_musculus.GRCm39.dna.primary_assembly.rdna.fa.gz"
    GENOMEFAtag=$(basename $GENOMEFA  .gz)
    GTFtag="Mus_musculus.GRCm39.104.rdna_rn18s.gtf"
    ;;
  GENCODEmm39|GENCODEm39)
    GENOMEDIR="$RPROF/data/mouse/gencode"
    GENOMEFA="GRCm39.primary_assembly.genome.rdna.fa.gz"
    GENOMEFAtag=$(basename $GENOMEFA  .gz)
    GTFtag="gencode.vM27.annotation.rdna_rn18s.gtf"
    ;;

  *)
    ;;
esac




GTFext=${GTFtag##*.}
GTFext=$(echo  "$GTFext"  |  tr  '[:upper:]' '[:lower:]' )

echo -e "GTFtag=[$GTFtag]"
echo -e "GTFext=[$GTFext]"



TMPPROJDIR="$TMPPROJDIR/$BASHPID"
TMPPROJDIRdata="$TMPPROJDIR/data"
TMPWDIR="$TMPPROJDIR/work"
TMPODIR="$TMPPROJDIR/out"
TMPPROJDIRstar="$TMPPROJDIR/STAR"
mkdir  -p  $TMPPROJDIR
mkdir  -p  $TMPPROJDIRdata
mkdir  -p  $TMPWDIR
mkdir  -p  $TMPODIR
mkdir  -p  $TMPPROJDIRstar


OUTTAG="STAR.index.$GENOMEFAtag--$GTFtag"

OUTDIR="$RPROF/data/STAR/$OUTTAG"
mkdir  -p  $OUTDIR

TMPOF="$TMPODIR/$OUTTAG"
TMPWF="$TMPWDIR/$OUTTAG"



cd $TMPPROJDIRdata







echo -e "copy forward...\n"

rsync  -tvh  $GENOMEDIR/$GENOMEFAtag.gz   $TMPPROJDIRdata
rsync  -tvh  $GENOMEDIR/$GTFtag.gz   $TMPPROJDIRdata

gunzip  $TMPPROJDIRdata/*.gz




echo -e "faidx...\n"

samtools  faidx  $TMPPROJDIRdata/$GENOMEFAtag






genomeLength=$(cat  $TMPPROJDIRdata/$GENOMEFAtag.fai  |  awk   'BEGIN {sum=0}  {sum+=$2}  END {print sum}' FS="\t"  OFS="\t")
echo -e "\t\t\tgenomeLength=[$genomeLength]\n"

STARflag_genomeSAindexNbases=""
if (( $genomeLength  <  2730871774 )); then

  # mouse genome m38: 2730871774   --genomeSAindexNbases 11
  echo -e "\tdetermine appropriate value for STAR flag --genomeSAindexNbases...\n"

  #2.2.5 Very small genome.
  #For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9, for 100 kiloBase genome, this is equal to 7.

  theCalc=$(perl -l -e "print int((log($genomeLength)/log(2))/2 - 1)")
  genomeSAflag=$(( theCalc<14?theCalc:14 ))
  #genomeSAflag=$((genomeSAflag - 3))
  genomeSAflag=$(( genomeSAflag<1?1:genomeSAflag ))
  # "Generally, --genomeSAindexNbases  needs to be scaled with the genome length, as ~min(14,log2(ReferenceLength)/2 - 1)"
  # https://groups.google.com/forum/#!msg/rna-star/9g8Uoe1Igho/b-NFBM02hhAJ

  #if [ "$genomeLength" -gt "10000000" ]; then
  #	genomeSAflag="9"
  #fi
  #if [ "$genomeLength" -gt "500000" ]; then
  #	genomeSAflag="7"
  #fi

  STARflag_genomeSAindexNbases="--genomeSAindexNbases $genomeSAflag" 
  echo -e "STARflag_genomeSAindexNbases=[$STARflag_genomeSAindexNbases]"

fi

readLength=100
nchr=$(wc -l < $TMPPROJDIRdata/$GENOMEFAtag.fai)
echo -e "\t\t\tNumberOfReferences=[$nchr]\n"

STARflag_genomeChrBinNbits=""
if (( $nchr  >  5000 )); then

  # mouse genome m38: nchr=66   --genomeChrBinNbits 18
  echo -e "\tdetermine appropriate value for STAR flag --genomeChrBinNbits...\n"

  # 2.2.6 Genome with a large number of references.
  #If you are using a genome with a large (>5,000) number of references (chrosomes/scaffolds), you may need to reduce the --genomeChrBinNbits to reduce RAM consumption. The following scaling is recommended: --genomeChrBinNbits = min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)]). For example, for 3 gigaBase genome with 100,000 chromosomes/scaffolds, this is equal to 15

  theCalc=$(Rscript -e "cat(sprintf('%.0f\n', min(18, log2(max($genomeLength/$nchr, $readLength)))))")
  genomeChrflag=$(( theCalc<1?1:theCalc ))

  STARflag_genomeChrBinNbits="--genomeChrBinNbits $genomeChrflag"
  echo -e "STARflag_genomeChrBinNbits=[$STARflag_genomeChrBinNbits]"

fi



STARgtfFlag=""
if [ "$GTFext" == "gff3"  ]; then
	STARgtfFlag="--sjdbGTFtagExonParentTranscript      Parent"
fi

echo -e "STARgtfFlag=[$STARgtfFlag]"




echo -e "\tgenerating STAR-format genome index file...\n"

STAR  \
    --runThreadN        $NUMTHREADS                         \
    --runMode           genomeGenerate                      \
    --genomeDir         $TMPPROJDIRstar                         \
    --genomeFastaFiles  $TMPPROJDIRdata/$GENOMEFAtag            \
    --outTmpDir         $TMPPROJDIR/star-tmp                    \
    --sjdbOverhang      $(($readLength-1))                  \
    --sjdbGTFfile       $TMPPROJDIRdata/$GTFtag                 \
    $STARgtfFlag  $STARflag_genomeSAindexNbases  $STARflag_genomeChrBinNbits



echo -e "\tcopying back...\n"
rsync  -av  $TMPPROJDIRstar/*  $OUTDIR


echo -e "\ttar...\b"

cd $RPROF/data/STAR/
tar  -zcv  -f $OUTTAG.tar.gz  -C $OUTTAG .



echo -e "\nDONE with  generate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh!\n\n"



