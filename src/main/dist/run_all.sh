# load gene models from NCBI for a number of species
. /etc/profile
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
cd $HOMEDIR

#              "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat"
SPECIES_LIST=( "1"     "2"     "3"   "6"   "5"      "7"        "4"          "9"   "13"     "14")

# old unreliable code: species common name can change
#SPECIES_LIST=( "human" "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat")


for SPECIES in "${SPECIES_LIST[@]}"; do
    $HOMEDIR/run_species.sh "$SPECIES"
done

# download gene_groups.gz file with gene-to-gene associations and load it into RGD_ASSOCIATIONS table
$HOMEDIR/load_gene_assoc.sh


