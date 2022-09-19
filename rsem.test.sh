#!/bin/bash 
#set  -eu  -o pipefail

source bashrc_rdna_variant_calling

date




#echo -e "\n\n---------------- rsem_eukaryotic_RNAseq_reads_to_annotated_genome.sh -----------------\n"
# Last updated:  2016-12-22

MIXEDFQ=""
INFQ1=""
INFQ2=""
LAYOUT=""
SPECIES=""
STRATEGY=""
PROJECT=""
INDIV=""
RUNID=""
KEEPRSEMTRANBAM=""
KEEPSTARBAM=""
STRANDED="reverse"
TAKEFQS=""
REPLACEFQS=""


OUTPREF=""
NUMTHREADS="1"
NUMTHREADS_for_cutadapt="1"
KEEPINTERM="0"

NODEMEM="8"

WORKFQDIR=`pwd`/work/temp-fastqs-var-call
mkdir -p $WORKFQDIR



while  test  $# -gt 0;  do
	case "$1"  in
		-h|--help|-help)
			echo -e "\n\nrsem_eukaryotic_RNAseq_reads_to_annotated_genome.sh\n"
			echo -e "This function will align RNA-seq reads to a genome using STAR and transcript annotations."
			echo -e "\n"

			echo -e "\nparameters:"
			echo -e "\t-u file\tinput single end fastq or fastq file of mixed mates for paired-end (will separate them appropriately). should be .gz"
			echo -e "\t-1 file\tinput fastq file for mate 1.  should be .gz"
			echo -e "\t-2 file\tinput fastq file for mate 2.  should be .gz"
			echo -e "\t-l layout\tlayout [paired, single]."
			echo -e "\t-s species\tthe species of the whole genome.  will load reference genome fasta and index if not provided (see \"-ix\" option)."
			echo -e "\t-strat  strategy\tseuqencing strategy, e.g. RNAseq or Riboseq."
			echo -e "\t-mem m\tamount of memory on node, in GB.  JAVA with flag -XmxAG uses A+3 G memory (there is a secret 3 Gb overhead).  Will tailor accordingly to node memory. [default: assumes m = 8]"

			echo -e "\n\tIf the following are provided, the output filename will be constructed and other parameters inferred."
			echo -e "\t-runid r\trunid"
			echo -e "\n"

			echo -e "\t-rsemTranBam 0/1\tkeep the transcriptome bam output from rsem [1] or discard [0]. Default: if strategy is ribosome profiling, then will keep, otherwise, discard."

			echo -e "\t-take\twill move the input fastq files to internal temporary storage."
			echo -e "\t-replace\twill replace fastq files to internal temporary storage, (if -take flag provided)."


			echo -e "\t-o file\toutput file prefix.  writes to <out>.bam"

			echo -e "\t-t pref\ttemporary file prefix."
			echo -e "\t-p\tnumber of threads for multithreading [default: 1]."
			echo -e "\t-k\tkeep intermediate files."

			echo -e "\nDESCRIPTION:"
			echo -e "\tThis function will align RNA-seq reads to the whole genome with STAR, accounting for transcript junctions, etc."

			echo -e "\nOUTPUT:"

			echo -e "\n\n"

			exit  0
			;;
		-u)
			# input file
			shift
			if test $# -gt 0; then
				MIXEDFQ="$1"
			else
				echo -e "\n\nERROR! mixed input fastq not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-1)
			# input mate 1 file
			shift
			if test $# -gt 0; then
				INFQ1="$1"
			else
				echo -e "\n\nERROR! input mate 1 fastq not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-2)
			# input mate 2 file
			shift
			if test $# -gt 0; then
				INFQ2="$1"
			else
				echo -e "\n\nERROR! input mate 2 fastq not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-l)
			# layout
			shift
			if test $# -gt 0; then
				LAYOUT="$1"
			else
				echo -e "\n\nERROR! layout not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-s)
			# species
			shift
			if test $# -gt 0; then
				SPECIES="$1"
			else
				echo -e "\n\nERROR! species not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-project)
			shift
			if test $# -gt 0; then
				PROJECT="$1"
			else
				echo -e "\n\nERROR! project not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-indiv)
			# indiv
			shift
			if test $# -gt 0; then
				INDIV="$1"
			else
				echo -e "\n\nERROR! indiv not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-runid)
			# runid
			shift
			if test $# -gt 0; then
				RUNID="$1"
			else
				echo -e "\n\nERROR! runid not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-mem)
			shift
			if test $# -gt 0; then
				NODEMEM="$1"
			else
				echo -e "\n\nERROR! vmem not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-rsemTranBam)
			KEEPRSEMTRANBAM="1"
			shift
			;;
                -starBam)
                        KEEPSTARBAM="1"
                        shift
                        ;;
		-take)
			TAKEFQS="-take"
			shift
			;;
		-replace)
			REPLACEFQS="1"
			shift
			;;
                -stranded)
                        shift
                        if test $# -gt 0; then
                                STRANDED="$1"
                        else
                                echo -e "\n\nERROR! stranded not specified.\n\n"
                                exit 1
                        fi
                        shift
                        ;;
		-o)
			# output prefix
			shift
			if test $# -gt 0; then
				OUTPREF="$1"
			else
				echo -e "\n\nERROR! output prefix not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-t)
			# temp
			shift
			if test $# -gt 0; then
				TMPWF="$1"
			else
				echo -e "\n\nERROR! temporary prefix not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-p)
			# threads
			shift
			if test $# -gt 0; then
				NUMTHREADS="$1"
				# for python3 environment
				NUMTHREADS_for_cutadapt="$1"
			else
				echo -e "\n\nERROR! number of threads not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-k)
			# keep intermediate files
			KEEPINTERM="1"
			shift
			;;
		*)
			echo -e "\n\nERROR!  nothing requested!\n"
			exit 1
			;;
	esac
done


if [ -z "$TMPPROJDIR"  -o  ! -d "$TMPPROJDIR" ]; then
	echo -e "\n\n\n\nERROR!  no temp dir!!!\n"
	exit 1
fi # tmpdir

TMPWF="$TMPPROJDIR/tmp.work-albqsr.$BASHPID"





############################ check parameters
echo -e "\n\n\n"
echo -e "MIXEDFQ=[$MIXEDFQ]"
echo -e "INFQ1=[$INFQ1]"
echo -e "INFQ2=[$INFQ2]"
echo -e "LAYOUT=[$LAYOUT]"
echo -e "SPECIES=[$SPECIES]"
echo -e "STRATEGY=[$STRATEGY]"
echo -e "PROJECT=[$PROJECT]"
echo -e "INDIV=[$INDIV]"
echo -e "RUNID=[$RUNID]"
echo -e "KEEPRSEMTRANBAM=[$KEEPRSEMTRANBAM]"
echo -e "STRANDED=[$STRANDED]"
echo -e "TAKEFQS=[$TAKEFQS]"
echo -e "REPLACEFQS=[$REPLACEFQS]"
echo -e "OUTPREF=[$OUTPREF]"
echo -e "NUMTHREADS=[$NUMTHREADS]"
echo -e "KEEPINTERM=[$KEEPINTERM]"
echo -e "NODEMEM=[$NODEMEM]"
echo -e "WORKFQDIR=[$WORKFQDIR]"
echo -e "\n\n\n"




TMPPROJDIRdata="$TMPWF.tmp-data-tmp"
mkdir  -p  $TMPPROJDIRdata

TMPPROJDIRfq="$TMPWF.tmp-fq"
mkdir  -p  $TMPPROJDIRfq





############### spack
#spack  load  -r  samtools@1.8
#spack  load  -r  r@3.5.0
#spack  load  -r  python@3.5.2



##################### error handling
### error hanbdling with trap
function replaceData {
        echo -e "\n\n\n"
        echo -e "error detected.  replacing data...\n"

	rsync  -tvh  $TMPPROJDIRfq/*   $WORKFQDIR

} # replaceData

if [ ! -z "$TAKEFQS" ]; then
	set  -E
	trap  'replaceData'  ERR  SIGINT  SIGTERM  SIGQUIT
fi # TAKEFQS





################# infer data

if [ ! -z "$RUNID"  -a  -z "$MIXEDFQ"  -a  -z "$INFQ1"  ]; then
	echo -e "runid is given.  inferring all data...\n"

	$RPROF/src/retrieve_subset_from_tabular_database.sh  -r $RUNID  -c project,individual,runid,species,layout,strategy  -o $TMPWF.info.sub

	tail  -n +2  $TMPWF.info.sub  >  $TMPWF.info.sub.tmp	    \
	&&  mv  $TMPWF.info.sub.tmp   $TMPWF.info.sub

	PROJECT=$(cat  $TMPWF.info.sub  |  cut -f 1)
	INDIV=$(cat  $TMPWF.info.sub  |  cut -f 2)
	SPECIES=$(cat  $TMPWF.info.sub  |  cut -f 4)
	LAYOUT=$(cat  $TMPWF.info.sub  |  cut -f 5)
	LAYOUT=$( echo "$LAYOUT"  |  tr '[:lower:]' '[:upper:]' )
	STRATEGY=$(cat  $TMPWF.info.sub  |  cut -f 6)

	if [ -z "$INFQ1" ]; then
		INFQ1="$TMPPROJDIRfq/$PROJECT.$INDIV.${RUNID}_1.fastq.gz"

		if [ "$LAYOUT" == "PAIRED"  -a  -z "$INFQ2" ]; then
			INFQ2="$TMPPROJDIRfq/$PROJECT.$INDIV.${RUNID}_2.fastq.gz"
		fi # LAYOUT
	fi # INFQ1

	if [ ! -s "$INFQ1" ]; then
		echo -e "getting raw fastq...\n"
                # begin of modification by H. Kim
                #$RPROF/src/get_raw_fastq_for_project.sh  $TAKEFQS  -runid $RUNID  -indiv $INDIV  -project $PROJECT   -o $TMPPROJDIRfq
                $RPROF/src/get_raw_fastq_for_project.sh  $TAKEFQS  -runid $RUNID  -indiv $INDIV  -project $PROJECT -l $LAYOUT  -o $TMPPROJDIRfq
                # end of modification
	fi # INFQ1

	ls -lrth $TMPPROJDIRfq

	if [ -z "$OUTPREF" ]; then
		OUTDIR="$RPROF/out/rsem-rnaseq-annotated-genome/$PROJECT/$INDIV/$RUNID"
		mkdir  -p  $OUTDIR
		OUTTAG="rsemRnaseq-annot.$PROJECT.$INDIV.$RUNID"
		OUTPREF="$OUTDIR/$OUTTAG"
	fi # OUTPREF

	echo -e "\n\n\n"
	echo -e "PROJECT=[$PROJECT]"
	echo -e "INDIV=[$INDIV]"
	echo -e "SPECIES=[$SPECIES]"
	echo -e "LAYOUT=[$LAYOUT]"
	echo -e "INFQ1=[$INFQ1]"
	echo -e "INFQ2=[$INFQ2]"
	echo -e "OUTDIR=[$OUTDIR]"
	echo -e "OUTTAG=[$OUTTAG]"
	echo -e "OUTPREF=[$OUTPREF]"
	echo -e "\n\n\n"

else  # proj, indiv
	OUTDIR=$OUTPREF

fi

LAYOUT=$( echo "$LAYOUT"  |  tr '[:lower:]' '[:upper:]' )

if [ "$LAYOUT" != "PAIRED"  -a  "$LAYOUT" != "SINGLE" ]; then
	echo -e "\n\n\n\nERROR!  bad layout given!!!\n"
	echo -e "LAYOUT=[$LAYOUT]\n\n"
	exit 1
fi


if [ -z "$KEEPRSEMTRANBAM"  -a  "$STRATEGY" == "Riboseq" ]; then
	echo -e "default: detected ribosome profiling, keeping rsem transcriptome bam...\n"
	KEEPRSEMTRANBAM="1"
fi # KEEPRSEMTRANBAM






if [ -z "$OUTDIR"  -o  ! -d "$OUTDIR" ]; then
	echo -e "\n\n\nERROR!   bad outdir given!!!\n"
	exit 1
fi #OUTDIR




TMPPROJDIRrsem="$TMPPROJDIR/star-ref"
mkdir  -p  $TMPPROJDIRrsem
TMPPROJDIRstar="$TMPPROJDIR/star-ref"
mkdir  -p  $TMPPROJDIRstar

OUTDIR=$(dirname  "$OUTPREF" )
OUTTAG=$(basename  "$OUTPREF" )
TMPODIR="$TMPWF.outdir"
mkdir  -p  $TMPODIR
TMPOF="$TMPODIR/$OUTTAG"







################################# copy forward
echo -e "copy forward...\n"

COPMIX=""
COPFQ1=""
COPFQ2=""

if [ ! -z "$MIXEDFQ" ]; then
	echo -e "\tmixed...\n"

	COPMIX=$( $RPROF/src/rsync_file_and_correct_suffix_appropriately.sh  -i $MIXEDFQ  -o $TMPWF.in.mixed.fq )
else
	echo -e "separate...\n"
	
	COPFQ1=$( $RPROF/src/rsync_file_and_correct_suffix_appropriately.sh  -i $INFQ1  -o $TMPWF.in.1.fq ) 

	if [ ! -z "$INFQ2" ]; then
		COPFQ2=$( $RPROF/src/rsync_file_and_correct_suffix_appropriately.sh  -i $INFQ2  -o $TMPWF.in.2.fq )
	fi
fi # copy





###################### prepare input fastq
if [ "$LAYOUT" == "PAIRED" ]; then
	if [ ! -z "$COPMIX" ]; then
		$RPROF/src/separate_mates_in_fastq.sh  -i $COPMIX  -o $TMPWF.sep  -t $TMPWF.tmp-sep-tmp  -p $NUMTHREADS
		COPFQ1="$TMPWF.sep.1.fq"
		COPFQ2="$TMPWF.sep.2.fq"
	fi
	
	if [ ! -s "$COPFQ1"  -o  ! -s "$COPFQ2" ]; then
		echo -e "\n\n\nERROR!  COPFQ1=[$COPFQ1]  or COPFQ2=[$COPFQ2] does not exist!!!\n"
		exit 1
	fi

fi # layout





if [ "$STRATEGY" == "Riboseq" ]; then
	echo -e "strategy is ribosome profiling, retrieving length limits..."

	if [ ! -z "$COPFQ2" ]; then
		echo -e "\n\n\nERROR!  ribosome profiling but copfq2 existds...\n"
		exit 1
	fi #COPFQ2
	## Start modification by Nitish
	$RPROF/src/retrieve_subset_from_tabular_database.sh  -d data2/database/161021.txt  -c project,individual,runid,species,layout,strategy,readLengthLower,readLengthUpper,adapter  -o $TMPWF.info.sub

	tail  -n +2  $TMPWF.info.sub  >  $TMPWF.info.sub.tmp	    \
	&&  mv  $TMPWF.info.sub.tmp   $TMPWF.info.sub

	READLENLOW=$(cat  $TMPWF.info.sub  |  cut -f 7)
	READLENUPP=$(cat  $TMPWF.info.sub  |  cut -f 8)
	ADAPTORSEQ=$(cat  $TMPWF.info.sub  |  cut -f 9)
	STRANDED="forward"

	echo -e "\n\n\n"
	echo -e "READLENLOW=[$READLENLOW]"
	echo -e "READLENUPP=[$READLENUPP]"
	echo -e "COPFQ1=[$COPFQ1]"
	echo -e "\n\n\n"


	### clip adapator
	echo -e "\tclipping and bounding...\n"

        # begin of modification
	#python3  $HOME/.local/bin/cutadapt   -a $ADAPTORSEQ  -m $READLENLOW  -M $READLENUPP  -o $TMPWF.trim.1.fq.gz   $COPFQ1  -j $NUMTHREADS
        # for 161021: $ADAPTORSEQ=CTGTAGGCACCATCAAT  $READLENLOW=18  $READLENUPP=42
	cutadapt   -a $ADAPTORSEQ  -m $READLENLOW  -M $READLENUPP  -o $TMPWF.trim.1.fq.gz   $COPFQ1  -j $NUMTHREADS_for_cutadapt
        # end of modification

	mv  $TMPWF.trim.1.fq.gz  $COPFQ1

else

  # begin of addition by H. Kim
  if [ "$LAYOUT" == "PAIRED" ]; then

    # Illumina TruSeq
    cutadapt   -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 18:18   -o $TMPWF.trim.1.fq.gz -p $TMPWF.trim.2.fq.gz  $COPFQ1 $COPFQ2  -j $NUMTHREADS_for_cutadapt
    mv -f $TMPWF.trim.1.fq.gz $COPFQ1
    mv -f $TMPWF.trim.2.fq.gz $COPFQ2
    
    #trim_galore --length 16 -q 20 --phred33 --paired --gzip $COPFQ1 $COPFQ2 -o $TMPPROJDIR
    #fname=$(echo "$COPFQ1" | sed -e "s/.fq.*//")
    #mv -f ${fname}_val_1.fq* $COPFQ1
    #fname=$(echo "$COPFQ2" | sed -e "s/.fq.*//")
    #mv -f ${fname}_val_2.fq* $COPFQ2

  else

    cutadapt   -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -m 18:18  -o $TMPWF.trim.1.fq.gz  $COPFQ1  -j $NUMTHREADS_for_cutadapt
    mv -f $TMPWF.trim.1.fq.gz $COPFQ1

    #trim_galore --length 16 -q 20 --phred33 --gzip $COPFQ1 -o $TMPPROJDIR
    #fname=$(echo "$COPFQ1" | sed -e "s/.fq.*//")
    #mv -f ${fname}_val_1.fq* $COPFQ1

  fi
  # end of addition

fi # STRATEGY





#################### load index

if [ -z "$SPECIES" ]; then
	echo -e "\n\n\nERROR!  SPECIES=[$SPECIES] specified!!!\n"
	exit 1
fi # SPECIES


## changes by Nitish
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





mkdir  -p  $TMPPROJDIRdata
rsync  -tvhL  $GENOMEDIR/$GENOMEFAtag.gz  $TMPPROJDIRdata
gunzip  $TMPPROJDIRdata/$GENOMEFAtag.gz

GTFbase=$(basename  $GTFtag)
rsync  -vhL  $GENOMEDIR/$GTFtag.gz                 $TMPPROJDIRdata
gunzip  -f $TMPPROJDIRdata/$GTFtag.gz

STARDIR="$RPROF/data/STAR"
STARTAG="STAR.index.$GENOMEFAtag--$GTFtag"

RSEMDIR="$RPROF/data/rsem-transcriptome"
RSEMTAG="$GENOMEFAtag--$GTFtag.rsem.transcriptome"





############ RSEM index
echo -e "rsem index...\n"
if [ ! -d $TMPPROJDIRrsem ] ; then
rsync  -tvh  $RSEMDIR/$RSEMTAG.tar.gz     $TMPPROJDIRrsem

echo -e "\tunpack index...\n"

tar  zxv  -C $TMPPROJDIRrsem  -f $TMPPROJDIRrsem/$RSEMTAG.tar.gz
fi



############ STAR index
echo -e "STAR index...\n"
if [ ! -d $TMPPROJDIRstar ]; then
rsync  -tvh  $STARDIR/$STARTAG.tar.gz   $TMPPROJDIRstar

echo -e "\tunpack index...\n"

tar  zxv  -C $TMPPROJDIRstar  -f $TMPPROJDIRstar/$STARTAG.tar.gz
fi




############################ RSEM
echo -e "calling rsem...\n"

PEFLAG=""
if [ "$LAYOUT" == "PAIRED" ]; then
	PEFLAG="--paired-end  "
fi # LAYOUT

STARPATH=$(dirname  $STAR)


BAMOUTFLAG="--no-bam-output"
if [ "$KEEPRSEMTRANBAM" == "1" ]; then
	BAMOUTFLAG=""
fi # KEEPRSEMTRANBAM

# begin of modification by H. Kim
if [ "$KEEPSTARBAM" == "1" ]; then
        BAMOUTFLAG="--output-genome-bam  \
                    --keep-intermediate-files"
fi
STRANDFLAG="--strandedness $STRANDED"
SEEDLENGTH=25
if [ "$STRATEGY" == "Riboseq" ]; then
  SEEDLENGTH=10
fi
# end of modification

echo -e "\n\n\n"
echo -e "PEFLAG=[$PEFLAG]"
echo -e "STARPATH=[$STARPATH]"
echo -e "BAMOUTFLAG=[$BAMOUTFLAG]"
# begin of modification by H. Kim
echo -e "STRANDFLAG=[$STRANDFLAG]"
echo -e "SEEDLENGTH=[$SEEDLENGTH]"
# end of modification
echo -e "\n\n\n"


rsem-calculate-expression			\
--star						\
--star-path	$STARPATH			\
--star-gzipped-read-file			\
--seed-length	$SEEDLENGTH			\
-p	$NUMTHREADS	 $BAMOUTFLAG	$STRANDFLAG	$PEFLAG	\
--temporary-folder	$TMPWF.rsem-tmp		\
--time						\
$COPFQ1   $COPFQ2				\
$TMPPROJDIRrsem/$RSEMTAG				\
$TMPOF


# begin of addition by H. Kim
if [ "$KEEPSTARBAM" == "1" ]; then
   echo -e "\t\tget mapped...\n"
   samtools  view  -hb  -F 4  $TMPOF.genome.bam  |  samtools  sort  -@ $NUMTHREADS  -O bam  -T $TMPWF.bam.tmp-tmp  -o $TMPOF.star_for_rsem.bam   -
   samtools  index  $TMPOF.star_for_rsem.bam
   mv -f  $TMPWF.rsem-tmp/${OUTTAG}Log.out $TMPOF.star_for_rsem.Log.out
   mv -f  $TMPWF.rsem-tmp/${OUTTAG}Log.progress.out $TMPOF.star_for_rsem.Log.progress.out
   mv -f  $TMPWF.rsem-tmp/${OUTTAG}Log.final.out $TMPOF.star_for_rsem.Log.final.out
else
   deleteIfExists   $TMPOF.star_for_rsem.bam
fi # KEEPSTARBAM
# end of addition


if [ "$KEEPRSEMTRANBAM" == "1" ]; then
	echo -e "\n\nsamtools sort transcript bam...\n"
	samtools  sort  -@ $NUMTHREADS  -T $TMPWF.sam-sort.tmp  -o $TMPOF.transcript.bam  $TMPOF.transcript.bam
	echo -e "index...\n"
	samtools  index  $TMPOF.transcript.bam
else
	deleteIfExists   $TMPOF.transcript.bam
fi # KEEPRSEMTRANBAM


echo -e "gzip...\n"
gzip  $TMPOF.*.results

find  $TMPODIR

echo -e "copy back...\n"
rsync  -tvh  $TMPOF.*   $OUTDIR

# begin of addition by H. Kim

#### STAR alignment

mkdir -p $TMPWF.star

# consider only uniquely mapped reads with --outFilterMultimapNmax 1
# default:
# --outFilterMultimapScoreRange default: 1
type_filter="unique"
type_sj="star_for_rsem"

ARG_FILTER=""
ARG_SJ=""
ARG_CHIM=""
ARG_PE=""

if [ "$STRATEGY" == "Riboseq" ]; then
   ARG_FILTER="--outFilterMultimapNmax              1 \
               --outFilterMismatchNmax            999 \
               --outFilterMismatchNoverReadLmax  0.04 \
               --outFilterMultimapScoreRange        1 \
	       --alignEndsType EndToEnd"
   # splicing switched off
   ARG_SJ="--alignIntronMax 1"

else

  case $type_filter in
    star_default)
       # STARmanual.pdf
       #     --outFilterMultimapNmax default: 10 int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as ”mapped to too many loci” in the Log.final.out .
       #    --outFilterMismatchNmax default: 10 int: alignment will be output only if it has no more mismatches than this value.
       #    --outFilterMismatchNoverReadLmax default: 1.0 real: alignment will be output onlyif its ratio of mismatches to *read* length is less than or equal to this value.
       #    --outFilterMultimapScoreRange default: 1 int: the score range below the maximum score for multimapping alignments
       ARG_FILTER="--outFilterMultimapNmax             10 \
                   --outFilterMismatchNmax             10 \
                   --outFilterMismatchNoverReadLmax   1.0 \
                   --outFilterMultimapScoreRange        1"
       ;;
    unique)
       #    --outFilterMismatchNoverReadLmax 0.04: max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.04*200=8 for the paired read
       ARG_FILTER="--outFilterMultimapNmax              1 \
                   --outFilterMismatchNmax            999 \
                   --outFilterMismatchNoverReadLmax  0.04 \
                   --outFilterMultimapScoreRange        1"
       ;;

    unique_endtoend)
       ARG_FILTER="--outFilterMultimapNmax              1 \
                   --outFilterMismatchNmax            999 \
                   --outFilterMismatchNoverReadLmax  0.04 \
                   --outFilterMultimapScoreRange        1 \
	           --alignEndsType EndToEnd"
       ;;

    encode)
       # STARmanual.pdf/ENCODE options
       ARG_FILTER="--outFilterMultimapNmax            20 \
                   --outFilterMismatchNmax           999 \
                   --outFilterMismatchNoverReadLmax 0.04 \
                   --outFilterMultimapScoreRange       1"
       ;;

    *)
       ARG_FILTER="--outFilterMultimapNmax           500 \
                   --outFilterMismatchNmax             2 \
                   --outFilterMismatchNoverReadLmax  1.0 \
                   --outFilterMultimapScoreRange       1"
       ;;
  esac

  case $type_sj in
    star_default)
        # STARmanual.pdf
        ARG_SJ="--outFilterType          Normal \
                --alignSJoverhangMin     5 \
                --alignIntronMin         21 \
                --alignSJDBoverhangMin 3 \
                --alignMatesGapMax 0 \
                --alignIntronMax 0 \
                --alignSJstitchMismatchNmax 0 -1 0 0"
        ;;
    encode)
        # STARmanual.pdf
        ARG_SJ="--outFilterType          BySJout \
                --alignSJoverhangMin     8 \
                --alignIntronMin         20 \
                --alignSJDBoverhangMin 1 \
                --alignMatesGapMax 1000000 \
                --alignIntronMax 1000000"
        ;;
    star-fusion)
        # https://github.com/STAR-Fusion/STAR-Fusion/wiki
        ARG_SJ="--alignSJDBoverhangMin 10 \
                --alignMatesGapMax 100000 \
                --alignIntronMax 100000 \
                --alignSJstitchMismatchNmax 5 -1 5 5"
        ;;
    star_for_rsem)
        ARG_SJ="--outFilterType          BySJout \
                --alignSJoverhangMin     8 \
                --alignIntronMin         20 \
                --alignSJDBoverhangMin 1 \
                --alignMatesGapMax 1000000 \
                --alignIntronMax 1000000"
        ;;
    *)
        # encode + star-fusion
        ARG_SJ="--outFilterType          BySJout \
                --alignSJoverhangMin     8 \
                --alignIntronMin         20 \
                --alignSJDBoverhangMin 10 \
                --alignMatesGapMax 100000 \
                --alignIntronMax 100000 \
                --alignSJstitchMismatchNmax 5 -1 5 5"
        ;;
  esac

fi # strategy


### Make changes in STAR path by Nitish
# STAR execution
STAR \
    --runThreadN        $NUMTHREADS \
    --runMode           alignReads  \
    --quantMode         TranscriptomeSAM GeneCounts \
    --genomeDir         $TMPPROJDIRstar  \
    --readFilesCommand  zcat \
    --readFilesIn       $COPFQ1 $COPFQ2 \
    --genomeLoad        NoSharedMemory \
    --sjdbGTFfile       $TMPPROJDIRdata/$GTFtag \
    $ARG_FILTER \
    --twopassMode       Basic \
    --outReadsUnmapped  Fastx \
    --outFileNamePrefix $TMPWF.star/ \
       --outSAMtype         BAM Unsorted \
       --outSAMunmapped     Within       \
       --outSAMstrandField  intronMotif  \
       --outSAMattrRGline   "ID:$PROJECT.$INDIV PL:illumina LB:RNAseq SM:$RUNID" \
    $ARG_SJ \
    $ARG_CHIM \
    $ARG_PE



OUTTAG="star-annot.$PROJECT.$INDIV.$RUNID"
TMPOF="$TMPODIR/$OUTTAG"

if [ "$KEEPSTARBAM" == "1" ]; then
   #cp -f $TMPWF.star/Aligned.out.bam $TMPOF.Aligned.out.bam
   #samtools  index  $TMPOF.Aligned.out.bam
   echo -e "\t\tget mapped...\n"
   samtools  view  -hb  -F 4  $TMPWF.star/Aligned.out.bam   | samtools  sort  -@ $NUMTHREADS  -O bam  -T $TMPWF.bam.tmp-tmp  -o $TMPOF.bam   -
   samtools  index  $TMPOF.bam
else
   #deleteIfExists   $TMPOF.Aligned.out.bam
   deleteIfExists   $TMPOF.bam
fi # KEEPSTARBAM


echo -e "\n\n\t\tget unmapped...\n"

#samtools view -f 4 $TMPWF.star/Aligned.out.bam > $TMPOF.unmapped.bam
if [ "$LAYOUT" == "PAIRED"  ]; then

  mv -f "$TMPWF.star/Unmapped.out.mate1"   "$TMPOF.unmapped_1.fastq"
  mv -f "$TMPWF.star/Unmapped.out.mate2"   "$TMPOF.unmapped_2.fastq"
  gzip -f $TMPOF.unmapped*.fastq

  #samtools collate -n 128 -u -O -@ $NUMTHREADS $TMPOF.unmapped.bam $TMPPROJDIR/bam-to-fastq-  |  samtools fastq -F 0x900 -c 6 -@ $NUMTHREADS -0 /dev/null -1 ${TMPOF}.unmapped_1.fastq.gz -2 ${TMPOF}.unmapped_2.fastq.gz -s ${TMPOF}.singleton.fastq.gz -

else

  mv -f "$TMPWF.star/Unmapped.out.mate1"   "$TMPOF.unmapped.fastq"
  gzip -f $TMPOF.unmapped.fastq

  #samtools collate -n 128 -u -O -@ $NUMTHREADS $TMPOF.unmapped.bam $TMPPROJDIR/bam-to-fastq-  |  samtools fastq -F 0x900 -c 6 -@ $NUMTHREADS -0 /dev/null -s ${TMPOF}.singleton.fastq.gz - > ${TMPOF}.unmapped.fastq.gz

fi


if [ -f $TMPWF.star/Chimeric.out.junction ]; then
  mv -f  $TMPWF.star/Chimeric.out.junction  $TMPOF.Chimeric.out.junction
fi
if [ -f $TMPWF.star/Chimeric.out.sam ]; then
  mv -f  $TMPWF.star/Chimeric.out.sam  $TMPOF.Chimeric.out.sam
fi
if [ -f $TMPWF.star/SJ.out.tab ]; then
  mv -f  $TMPWF.star/SJ.out.tab  $TMPOF.SJ.out.tab
fi
mv -f  $TMPWF.star/Aligned.toTranscriptome.out.bam $TMPOF.Aligned.toTranscriptome.out.bam
mv -f  $TMPWF.star/ReadsPerGene.out.tab $TMPOF.ReadsPerGene.out.tab
mv -f  $TMPWF.star/Log.out $TMPOF.Log.out
mv -f  $TMPWF.star/Log.progress.out $TMPOF.Log.progress.out
mv -f  $TMPWF.star/Log.final.out $TMPOF.Log.final.out


echo -e "copy back...\n"
rsync  -tvh  $TMPOF.*   $OUTDIR

# make links
ln  -sf  $OUTDIR/$OUTTAG.bam  $OUTDIR/$RUNID.bam
ln  -sf  $OUTDIR/$OUTTAG.bam.bai  $OUTDIR/$RUNID.bam.bai

# bigwig
bamCoverage --normalizeUsing None -b $OUTDIR/$RUNID.bam -o $OUTDIR/$RUNID.bw


### HTSEQ
echo -e "\nhtseq-count...\n"

OUTTAG="htseq-annot.$PROJECT.$INDIV.$RUNID"
TMPOF="$TMPODIR/$OUTTAG"

if [ $STRANDED == "none" ]; then STRANDED="no"; fi
if [ $STRANDED == "forward" ]; then STRANDED="yes"; fi

STRANDFLAG="--stranded $STRANDED"
echo -e "STRANDFLAG=[$STRANDFLAG]"

htseq-count --quiet -r name  $STRANDFLAG  -f bam -m union -t exon $TMPWF.star/Aligned.out.bam $TMPPROJDIRdata/$GTFtag -a 10 -i gene_id > $TMPOF.genes.results




echo -e "gzip...\n"
gzip  $TMPOF.*.results

find  $TMPODIR

echo -e "copy back...\n"
rsync  -tvh  $TMPOF.*   $OUTDIR


echo -e "\n\nDONE!\n\n\n\n"

date
