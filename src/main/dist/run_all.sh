# script to run first rat EntrezGene pipeline followed by Ortholog loading and Ortholog ftp export
# and then human and mouse EntrezGene pipelines
. /etc/profile
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
cd $HOMEDIR

SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" )

for SPECIES in "${SPECIES_LIST[@]}"; do
    $HOMEDIR/run_species.sh "$SPECIES"
done

# download gene_groups.gz file with gene-to-gene associations and load it into RGD_ASSOCIATIONS table
$HOMEDIR/load_gene_assoc.sh

# download gene_history.gz file and withdraw or merge genes for all species except rat
$HOMEDIR/handle_ncbi_gene_history.sh

