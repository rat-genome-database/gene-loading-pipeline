#download and process all the genes for given species and map key
# which do not have strand information available
# parameters: 'species' 'map_key'
# f.e. no_strand.sh rat 360
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu

echo  "starting $1 EntrezGene pipeline in no_strand mode"
cd $HOMEDIR
java -Dspring.config=../properties/default_db.xml \
    -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
    -jar lib/EntrezGeneLoading.jar \
    -no_strand $2 \
    -species $1 > no_strand_$1_$2.log
mailx -s "[$SERVER] $1 EntrezGene pipeline finished running" $ELIST < no_strand_$1_$2.log