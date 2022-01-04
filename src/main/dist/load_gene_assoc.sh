#load gene-to-gene associations from gene_group.gz file
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
LOGFILE=$HOMEDIR/assoc.log

ELIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
    ELIST="rgd.devops@mcw.edu,rgd.pipelines@mcw.edu"
fi

cd $HOMEDIR
echo "" > $LOGFILE

#              "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat"
SPECIES_LIST=( "1"     "2"     "3"   "6"   "5"      "7"        "4"          "9"   "13"     "14")

# old unreliable code: species common name can change
#SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat")

for SPECIES in "${SPECIES_LIST[@]}"; do

    echo  "load gene-to-gene associations for $SPECIES"
    java -Dspring.config=../properties/default_db2.xml \
        -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
        -jar lib/EntrezGeneLoading.jar \
        -load_gene_associations \
        -species "$SPECIES" >> $LOGFILE
done

echo "load gene associations complete"
mailx -s "[$SERVER] EntrezGene pipeline gene associations loaded" $ELIST < $LOGFILE
