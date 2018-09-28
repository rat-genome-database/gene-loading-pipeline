# script to run first rat EntrezGene pipeline followed by Ortholog loading and Ortholog ftp export
# and then human and mouse EntrezGene pipelines
. /etc/profile
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
cd $HOMEDIR

$HOMEDIR/run_species.sh rat
$HOMEDIR/run_species.sh mouse
$HOMEDIR/run_species.sh human
$HOMEDIR/run_species.sh chinchilla
$HOMEDIR/run_species.sh bonobo
$HOMEDIR/run_species.sh dog
$HOMEDIR/run_species.sh squirrel

# download gene_groups.gz file with gene-to-gene associations and load it into RGD_ASSOCIATIONS table
$HOMEDIR/load_gene_assoc.sh

# download gene_history.gz file and withdraw or merge genes for all species except rat
$HOMEDIR/handle_ncbi_gene_history.sh

