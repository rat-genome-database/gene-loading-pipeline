#on demand load of all mitochondrial genes
#
echo  "starting EntrezGene pipeline in mitochondrial mode"
cd /home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" )

for SPECIES in "${SPECIES_LIST[@]}"; do
    $HOMEDIR/run_species.sh -mitochondrial "$SPECIES" > mt${SPECIES}.log
    mailx -s "[$SERVER] $SPECIES EntrezGene pipeline for mitochondrial genes OK" mtutaj@mcw.edu < mt${SPECIES}.log
done

