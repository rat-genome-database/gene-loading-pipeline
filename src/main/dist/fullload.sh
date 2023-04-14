#for given species, download and process all genes that have been modified since 2000/01/01 at NCBI
# parameter is common species name: rat|mouse|human|...
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
TODAY=`date +%Y-%m-%d`
LOGFILE=fullload_$1_${TODAY}.log

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="mtutaj@mcw.edu"
fi

echo  "starting full-load Gene pipeline for $1"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
    -jar lib/EntrezGeneLoading.jar \
    -download+process fullload auto \
    -species $1 > $LOGFILE
mailx -s "[$SERVER] Full Load Gene pipeline finished running for $1" $ELIST < $LOGFILE


