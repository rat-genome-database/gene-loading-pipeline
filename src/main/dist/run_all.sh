# load gene models from NCBI for a number of species
. /etc/profile
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
cd $HOMEDIR

SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat" )

for SPECIES in "${SPECIES_LIST[@]}"; do
    $HOMEDIR/run_species.sh "$SPECIES"
done

# download gene_groups.gz file with gene-to-gene associations and load it into RGD_ASSOCIATIONS table
$HOMEDIR/load_gene_assoc.sh


