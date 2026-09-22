#runs a load for all active biological regions of given species
# (objects of type 'biological region', and genes still typed 'biological-region', which the load converts)
# note: run_all_genes.sh covers genes only, biological regions are not part of it
# parameter is rat|mouse|human
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="$ELIST rgd.pipelines@mcw.edu"
fi

echo  "starting $1 EntrezGene pipeline for biological regions"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
    -jar lib/EntrezGeneLoading.jar \
    -all_biological_regions \
    -species "$1" > "$1_all_biological_regions.log"
mailx -s "[$SERVER] $1 EntrezGene pipeline for biological regions finished running" $ELIST < "$1_all_biological_regions.log"
