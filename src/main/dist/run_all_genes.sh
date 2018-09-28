#runs a load for all active genes of given species
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
    -all_genes \
    -species $1 > $1_all_genes.log
mailx -s "[$SERVER] $1 EntrezGene pipeline finished running" $ELIST < $1_all_genes.log

#mail a file with changed gene symbols to Stan
mailx -s "[$SERVER] $1 gene symbol conflicts" mtutaj@mcw.edu,slaulederkind@mcw.edu < logs/symbols.log

