#!/bin/bash

echo -e "\n\n---------------- generate_rsem_transcriptome_index.sh -----------------\n"
echo -e "This script will generate a reference index for RSEM.\n"
echo -e "\t\t command :: generate_rsem_transcriptome_index.sh mm39\n"


GENOMETAG="$1"

echo -e "\t\t\tGENOMETAG=[$GENOMETAG]\n"

SPECIES=$GENOMETAG
### Modification by Nitish
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

RSEMDIR="$RPROF/data/rsem-transcriptome"
RSEMTAG="$GENOMEFAtag--$GTFtag.rsem.transcriptome"

mkdir  -p  $RSEMDIR

TMPPROJDIR="$TMPPROJDIR/$BASHPID"
TMPPROJDIRdata="$TMPPROJDIR/data"
TMPODIR="$TMPPROJDIR/out"
TMPWDIR="$TMPPROJDIR/work"
mkdir  -p  $TMPPROJDIR
mkdir  -p  $TMPPROJDIRdata
mkdir  -p  $TMPODIR
mkdir  -p  $TMPWDIR


OUTTAG="$GENOMEFAtag--$GTFtag.rsem.transcriptome"
# end of modification
TMPOF="$TMPODIR/$OUTTAG"

################################################## parameters
RSEMNGVECK="15"


cd  $TMPPROJDIRdata

echo -e "\tcopy forward...\n"

rsync  -vhL  $GENOMEDIR/$GTFtag.gz  $TMPPROJDIRdata
rsync  -tvh  $GENOMEDIR/$GENOMEFAtag.gz   $TMPPROJDIRdata
# end of modification
gunzip  $TMPPROJDIRdata/*.gz


echo -e "\tgenerating new rsem reference...\n"

# begin of modification by Nitish
rsem-prepare-reference          \
--gtf  $TMPPROJDIRdata/$GTFtag            \
--num-threads 10                          \
$TMPPROJDIRdata/$GENOMEFAtag              \
$TMPOF
# end of modification



echo -e "\tmake ngvector...\n"

$RSEM/rsem-generate-ngvector					\
-k				$RSEMNGVECK			\
$TMPOF.transcripts.fa						\
$TMPOF



echo -e "\tcompress...\n"

cd  $TMPODIR
tar  -zcvf  $OUTTAG.tar.gz   $OUTTAG.*
cd  $TMPPROJDIRData


echo -e "\tcopying back...\n"
rsync  -tvh  $TMPOF.tar.gz   $RSEMDIR



echo -e "\nDONE with  generate_rsem_transcriptome_index.sh!\n\n"


