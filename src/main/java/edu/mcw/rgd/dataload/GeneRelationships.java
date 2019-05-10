package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since 10/26/11
 * The contents of gene_group.gz file (as explained in ftp://ftp.ncbi.nih.gov/gene/README)
 * <pre>
   gene_group - recalculated daily
   report of genes and their relationships to other genes

   tab-delimited
   one line per GeneID
   Column header line is the first line in the file.

   NOTE: This file is not comprehensive, and contains
   a subset of information summarizing gene-gene relationships.

   Please consider HomoloGene and ProteinClusters
   as additional sources of information.

   ftp://ftp.ncbi.nih.gov/pub/HomoloGene/
   ftp://ftp.ncbi.nih.gov/genomes/Bacteria/CLUSTERS/

   Relationships are reported symmetrically, and currently include:
       Potential readthrough sibling
       Readthrough child
       Readthrough parent
       Readthrough sibling
       Related functional gene
       Related pseudogene
       Region member
       Region parent

    'Ortholog' relationships have been discontinued by NCBI as of beginning 2018
    (Ortholog relationships are now available from file gene_orthologs.gz, which is processed by Ortholog pipeline)
 </pre>
 */
public class GeneRelationships {

    private String version;
    private String geneGroupFile;

    // list of all associations for assoc types: readthrough_gene | related_functional_gene | related_pseudogene;
    // during loading, if incoming association is matching RGD, its assoc key is removed from this set;
    // when processing is done, assocs that have been left will be removed from db entirely
    final private Map<Integer, Association> assocKeys = new HashMap<>();

    final private Set<Integer> supportedSpeciesTaxIds = new HashSet<>();

    public void loadAssociations(final int speciesTypeKey) throws Exception {

        System.out.println(getVersion());

        final EGDAO dao = EGDAO.getInstance();

        // get tax id from species type key
        final int taxId = SpeciesType.getTaxonomicId(speciesTypeKey);

        // load tax ids for supported type keys
        for( int supportedSpeciesTypeKey: SpeciesType.getSpeciesTypeKeys() ) {
            supportedSpeciesTaxIds.add(SpeciesType.getTaxonomicId(supportedSpeciesTypeKey));
        }

        // download the file
        final String localFile = this.downloadAssocFile();

        // read all assoc keys for assoc with types readthrough_gene | related_functional_gene | related_pseudogene
        readAssocKeys(speciesTypeKey);

        // parses all of the records from gene_group.gz file and turns them into a list of Association objects
        List<AssocRecord> incomingAssocs = new ArrayList<>();
        BufferedReader reader = Utils.openReader(localFile);
        String line;
        while( (line=reader.readLine()) != null ) {
            // skip comments
            if( line.startsWith("#") )
                continue;
            String[] cols = line.split("\\t", -1);

            AssocRecord rec = new AssocRecord();
            rec.taxId1 = Integer.parseInt(cols[0]);
            rec.geneId1 = Integer.parseInt(cols[1]);
            rec.relType = cols[2];
            rec.taxId2 = Integer.parseInt(cols[3]);
            rec.geneId2 = Integer.parseInt(cols[4]);

            incomingAssocs.add(rec);
        }
        reader.close();

        // to minimize conflicts, shuffle the incomning data
        Collections.shuffle(incomingAssocs);

        CounterPool counters = new CounterPool();
        incomingAssocs.parallelStream().forEach( rec -> {

            // matches eg ids against rgd
            // if association species for master rgd id is different from the current species, skip it
            // the second species must be among the supported species
            if( rec.taxId1!=taxId || !supportedSpeciesTaxIds.contains(rec.taxId2) ) {
                counters.increment("LOAD_DIFF_SPECIES"); // different species: skip this record during loading
                return;
            }

            if( rec.taxId1==taxId && rec.taxId2!=taxId
                    || rec.taxId1!=taxId && rec.taxId2==taxId ) {
                counters.increment("ORTHO_ASSOC_FOUND");
            }

            if( rec.taxId1==rec.taxId2 ) {
                counters.increment("SPECIES_SAME");
            } else {
                counters.increment("SPECIES_DIFF");
            }

            try {
                rec.rgdId1 = dao.matchRgdIdByEgId(rec.geneId1);
                rec.rgdId2 = dao.matchRgdIdByEgId(rec.geneId2);

                if( !(rec.rgdId1>0 && rec.rgdId2>0) ) {
                    counters.increment("LOAD_SKIP (no match against RGD)");
                    return;
                }

                Association assoc = new Association();
                switch (rec.relType) {
                    case "Potential readthrough sibling":
                        assoc.setAssocType("readthrough_gene");
                        assoc.setAssocSubType("potential_readthrough_sibling");
                        break;
                    case "Readthrough child":
                        assoc.setAssocType("readthrough_gene");
                        assoc.setAssocSubType("readthrough_child");
                        break;
                    case "Readthrough parent":
                        assoc.setAssocType("readthrough_gene");
                        assoc.setAssocSubType("readthrough_parent");
                        break;
                    case "Readthrough sibling":
                        assoc.setAssocType("readthrough_gene");
                        assoc.setAssocSubType("readthrough_sibling");
                        break;
                    case "Related functional gene":
                        assoc.setAssocType("related_functional_gene");
                        break;
                    case "Related pseudogene":
                        assoc.setAssocType("related_pseudogene");
                        break;
                    case "Region member":
                        assoc.setAssocType("region");
                        assoc.setAssocSubType("region_member");
                        break;
                    case "Region parent":
                        assoc.setAssocType("region");
                        assoc.setAssocSubType("region_parent");
                        break;
                    case "Ortholog":
                        counters.increment("LOAD_ORTHOLOG_SKIP"); // skip this record during loading
                        return;
                    default:
                        counters.increment("LOAD_UNEXPECTED_GENE_RELATION_TYPE_" + rec.relType);
                        return;
                }

                assoc.setMasterRgdId(rec.rgdId1);
                assoc.setDetailRgdId(rec.rgdId2);
                assoc.setSrcPipeline("ENTREZGENE");
                rec.assoc = assoc;

                // now determine the action for the loader:
                // LOAD_MATCH -- perfect match, nothing to do
                // LOAD_INSERT -- new association to be inserted
                List<Association> assocList = dao.getAssociationsByRgdId(assoc.getMasterRgdId());
                if( assocList.contains(assoc) ) {

                    int matchingAssocIndex = assocList.indexOf(assoc);
                    Association matchingAssocInRgd = assocList.get(matchingAssocIndex);

                    // now check if subtype did not change (subtype is not part of an unique index)
                    if( matchingAssocInRgd.equalsWithSubType(assoc) ) {
                        counters.increment("LOAD_MATCH");
                    }
                    else {
                        // the assoc has to be updated -- subtype has been changed
                        assoc.setAssocKey(matchingAssocInRgd.getAssocKey());
                        assoc.setCreationDate(matchingAssocInRgd.getCreationDate());
                        counters.increment("LOAD_UPDATE");

                        dao.updateAssociation(rec.assoc);
                    }

                    synchronized(assocKeys) {
                        assocKeys.remove(matchingAssocInRgd.getAssocKey());
                    }
                }
                else {
                    counters.increment("LOAD_INSERT");

                    dao.insertAssociation(rec.assoc);
                }
            } catch(Exception e) {
                throw new RuntimeException(e);
            }
        });

        // report associations not handled by the pipeline
        int deletedAssocCount = 0;
        for( Association assoc: assocKeys.values() ) {
            deletedAssocCount += dao.deleteAssociation(assoc);
        }
        counters.add("LOAD_DELETE", deletedAssocCount);

        System.out.println(counters.dumpAlphabetically());

        System.out.println("OK!");
    }

    // download file with associations; return name of local copy of this file
    private String downloadAssocFile() throws Exception {

        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(getGeneGroupFile());
        downloader.setLocalFile("data/gene_group.gz");
        downloader.setPrependDateStamp(true);
        downloader.setUseCompression(true);
        return downloader.download();
    }

    private void readAssocKeys(int speciesTypeKey) throws Exception {

        EGDAO dao = EGDAO.getInstance();
        readAssocKeys(dao, "readthrough_gene", speciesTypeKey);
        readAssocKeys(dao, "related_functional_gene", speciesTypeKey);
        readAssocKeys(dao, "related_pseudogene", speciesTypeKey);
        readAssocKeys(dao, "ortholog", speciesTypeKey);
        readAssocKeys(dao, "region", speciesTypeKey);
    }

    private void readAssocKeys(EGDAO dao, String assocType, int speciesTypeKey) throws Exception {

        for(Association assoc: dao.getAssociationsByType(assocType, speciesTypeKey) ) {
            this.assocKeys.put(assoc.getAssocKey(), assoc);
        }
    }

    public String getGeneGroupFile() {
        return geneGroupFile;
    }

    public void setGeneGroupFile(String geneGroupFile) {
        this.geneGroupFile = geneGroupFile;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public class AssocRecord {

        // incoming data
        int taxId1;
        int geneId1;
        String relType; // relationship type
        int taxId2;
        int geneId2;

        // computed data
        int rgdId1;
        int rgdId2;
        Association assoc; // null if association cannot be created
    }
}
