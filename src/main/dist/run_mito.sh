#on demand load of all mitochondrial genes
#
echo  "starting EntrezGene pipeline in mitochondrial mode"
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
cd $HOMEDIR
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" )

for SPECIES in "${SPECIES_LIST[@]}"; do
    java -Dspring.config=../properties/default_db2.xml \
        -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
        -jar lib/EntrezGeneLoading.jar -mitochondrial -species "$SPECIES" > mt${SPECIES}.log
    mailx -s "[$SERVER] mitochondrial $SPECIES EntrezGene pipeline OK" mtutaj@mcw.edu < mt${SPECIES}.log
done

