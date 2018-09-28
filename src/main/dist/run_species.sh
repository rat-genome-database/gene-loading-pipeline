#runs weekly loading for one species
# parameter is rat|mouse|human
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=rgd.developers@mcw.edu
if [ "$SERVER" == "REED" ]
then
    ELIST="$ELIST,rgd.pipelines@mcw.edu"
fi

echo  "starting $1 EntrezGene pipeline"
cd $HOMEDIR
java -Dspring.config=../properties/default_db.xml \
    -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
    -jar EntrezGeneLoading.jar \
    -download+process auto auto \
    -species $1 > $1.log
mailx -s "[$SERVER] $1 EntrezGene pipeline finished running" $ELIST < $1.log


