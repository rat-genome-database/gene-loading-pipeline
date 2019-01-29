#download gene_history.gz file from NCBI containing withdrawn and secondary NCBI gene ids
# and process it for all species except rat: withdraw genes or merge them (for secondary NCBI gene ids)
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
YMD=`date +"%Y-%m-%d"`
LOGFILE="$HOMEDIR/${YMD}_ncbi_gene_history.log"

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="mtutaj@mcw.edu,jrsmith@mcw.edu,slaulederkind@mcw.edu"
fi

cd $HOMEDIR

echo "handle ncbi gene history"
java -Dspring.config=../properties/default_db.xml \
    -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
    -jar lib/EntrezGeneLoading.jar \
    -ncbi_gene_history 2>&1 > $LOGFILE

echo "ncbi gene history file processed"
mailx -s "[$SERVER] NcbiGene pipeline gene history complete" $ELIST < $LOGFILE
