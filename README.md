# gene-loading-pipeline

Loads gene models from NCBI (genes and transcripts) for all species in RGD
 (rat, mouse, human, chinchilla, bonobo, dog, pig, squirrel, green monkey, naked mole rat, ...)

Created genes have their gene source set to 'NCBI'.

There is a separate logic for rat genes, and for non-rat genes.

Matching logic for non-rat genes (to determine if an incoming gene matches a gene in RGD):

    1. match by NCBI Gene Id (there could be only one active gene with a given NCBI Gene id)
    2. if not found, match by HGNC/MGI id -- only for human/mouse genes
    3. if not found, match by Ensembl Gene Id
    
## Nomenclature changes for non-rat genes.
    
    1. If the gene nomenclature source authority is 'NCBI', if there are any changes to gene name or symbol,
    then gene name and/or symbol are updated and nomenclature event is generated
    2. If the gene nomenclature source authority is not 'NCBI', then any changes
    to gene name and/or symbol are suppressed       