#!/bin/bash
# loads, or restores, the transcripts of one assembly from an NCBI GFF3 file,
# f.e. from an archived annotation release of an older assembly:
#
#   load_transcripts_from_gff3.sh <map_key> <gff3_file>
#   load_transcripts_from_gff3.sh 372 /data/GCF_015227675.2_mRatBN7.2_genomic.gff.gz
#
# transcripts already in RGD are matched by accession; transcripts detached in the past are restored
# under their old rgd id (per STABLE_TRANSCRIPTS); existing feature objects are bound, not duplicated;
# genes on unplaced scaffolds and genes inactive in RGD are skipped
#
if [ $# -ne 2 ]; then
    echo "usage: $0 <map_key> <gff3_file>"
    exit 1
fi
MAP_KEY="$1"
GFF3_FILE="$2"

if ! [[ "$MAP_KEY" =~ ^[0-9]+$ ]]; then
    echo "map_key must be a number: $MAP_KEY"
    exit 1
fi
if [ ! -f "$GFF3_FILE" ]; then
    echo "gff3 file not found: $GFF3_FILE"
    exit 1
fi
# absolute path: the loader runs from HOMEDIR
GFF3_FILE=$(readlink -f "$GFF3_FILE")

HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="$ELIST rgd.pipelines@mcw.edu"
fi

LOG="transcripts_gff3_map$MAP_KEY.log"

echo "starting transcript gff3 loader: map_key=$MAP_KEY file=$GFF3_FILE"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
    -jar lib/EntrezGeneLoading.jar \
    -transcripts_from_gff3 "$MAP_KEY" "$GFF3_FILE" \
    -species rat > "$LOG" 2>&1

# the log lists every gene processed; the final counters are at its end
tail -100 "$LOG" | mailx -s "[$SERVER] transcript gff3 loader for map_key $MAP_KEY finished running" $ELIST
