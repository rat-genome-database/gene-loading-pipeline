#runs weekly loading for one species
# parameter is rat|mouse|human|...
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="mtutaj@mcw.edu"
fi

echo  "starting full-load Gene pipeline for %1"
cd $HOMEDIR
java -Dspring.config=../properties/default_db.xml \
    -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
    -jar lib/EntrezGeneLoading.jar \
    -download+process 2000/01/01 auto \
    -species $1 > fullload_$1.log
mailx -s "[$SERVER] Full Load Gene pipeline finished running for %1" $ELIST < fullload_$1.log


