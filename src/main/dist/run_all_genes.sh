#runs a load for all active genes of given species
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
    -all_genes \
    -species $1 > $1_all_genes.log
mailx -s "[$SERVER] $1 EntrezGene pipeline finished running" $ELIST < $1_all_genes.log

#mail a file with changed gene symbols to Stan
mailx -s "[$SERVER] $1 gene symbol conflicts" mtutaj@mcw.edu,slaulederkind@mcw.edu < logs/symbols.log

