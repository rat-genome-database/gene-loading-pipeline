#runs weekly loading for one species
# parameter is rat|mouse|human
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="$ELIST,rgd.pipelines@mcw.edu"
fi

echo  "starting $1 EntrezGene pipeline"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
    -jar lib/EntrezGeneLoading.jar \
    -download+process auto auto \
    -species $1 > $1.log
mailx -s "[$SERVER] $1 EntrezGene pipeline finished running" $ELIST < $1.log


