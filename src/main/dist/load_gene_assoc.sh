#load gene-to-gene associations from gene_group.gz file
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
LOGFILE=$HOMEDIR/assoc.log

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="rgd.developers@mcw.edu,rgd.pipelines@mcw.edu"
fi

cd $HOMEDIR
echo "" > $LOGFILE

SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" )

for SPECIES in "${SPECIES_LIST[@]}"; do

    echo  "load gene-to-gene associations for $SPECIES"
    java -Dspring.config=../properties/default_db.xml \
        -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
        -jar EntrezGeneLoading.jar \
        -load_gene_associations \
        -species "$SPECIES" >> $LOGFILE
done

echo "load gene associations complete"
mailx -s "[$SERVER] EntrezGene pipeline gene associations loaded" $ELIST < $LOGFILE
