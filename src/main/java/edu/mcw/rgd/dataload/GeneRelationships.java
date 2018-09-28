package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.pipelines.PipelineManager;
import edu.mcw.rgd.pipelines.PipelineRecord;
import edu.mcw.rgd.pipelines.RecordPreprocessor;
import edu.mcw.rgd.pipelines.RecordProcessor;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.PipelineLogger;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.Map;
import java.util.zip.GZIPInputStream;

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

    PipelineLogger dbLog = PipelineLogger.getInstance();

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

        // first processor parses all of the records from gene_group.gz file
        // and turns them into a bunch of Association objects
        PipelineManager manager = new PipelineManager();

        manager.addPipelineWorkgroup(new RecordPreprocessor() {
            @Override
            public void process() throws Exception {

                int recNo = 0;
                BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(localFile))));
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

                    rec.setRecNo(++recNo);
                    getSession().putRecordToFirstQueue(rec);
                }
                reader.close();
            }
        }, "PP", 1, 0);

        // 2nd processor: matches eg ids against rgd
        manager.addPipelineWorkgroup(new RecordProcessor() {
            @Override
            public void process(PipelineRecord pipelineRecord) throws Exception {

                AssocRecord rec = (AssocRecord) pipelineRecord;

                rec.getFlags().clear();

                // if association species for master rgd id is different from the current species, skip it
                if( rec.taxId1!=taxId ) {
                    rec.setFlag("LOAD_DIFF_SPECIES"); // different species: skip this record during loading
                    return;
                }
                // the second species must be among the supported species
                if( !supportedSpeciesTaxIds.contains(rec.taxId2) ) {
                    rec.setFlag("LOAD_DIFF_SPECIES"); // different species: skip this record during loading
                    return;
                }

                if( rec.taxId1==taxId && rec.taxId2!=taxId
                 || rec.taxId1!=taxId && rec.taxId2==taxId ) {
                    rec.setFlag("ORTHO_ASSOC_FOUND");
                }

                if( rec.taxId1==rec.taxId2 ) {
                    rec.setFlag("SPECIES_SAME");
                } else {
                    rec.setFlag("SPECIES_DIFF");
                }

                rec.rgdId1 = dao.matchRgdIdByEgId(rec.geneId1);
                rec.rgdId2 = dao.matchRgdIdByEgId(rec.geneId2);

                // create a new assoc object
                if( rec.rgdId1>0 && rec.rgdId2>0 ) {
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
                            rec.setFlag("LOAD_ORTHOLOG_SKIP"); // skip this record during loading
                            return;
                        default:
                            System.out.println("unknown gene relation type: " + rec.relType);
                            rec.setFlag("LOAD_UNEXPECTED");
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
                            rec.setFlag("LOAD_MATCH");
                        }
                        else {
                            // the assoc has to be updated -- subtype has been changed
                            assoc.setAssocKey(matchingAssocInRgd.getAssocKey());
                            assoc.setCreationDate(matchingAssocInRgd.getCreationDate());
                            rec.setFlag("LOAD_UPDATE");
                        }
                        assocKeys.remove(matchingAssocInRgd.getAssocKey());
                    }
                    else {
                        rec.setFlag("LOAD_INSERT");
                    }
                }
                else {
                    //System.out.println("skipping: no match EGID:"+rec.geneId1+"<==>RGDID:"+rec.rgdId1+", EGID:"+rec.geneId2+"<==>RGDID:"+rec.rgdId2);
                    rec.setFlag("LOAD_SKIP"); // skip this record during loading
                }
            }
        }, "QC", 5, 0);

        manager.addPipelineWorkgroup(new RecordProcessor() {
            @Override
            public void process(PipelineRecord pipelineRecord) throws Exception {

                AssocRecord rec = (AssocRecord) pipelineRecord;

                if( rec.isFlagSet("ORTHO_ASSOC_FOUND") ) {
                    getSession().incrementCounter("ORTHO_ASSOC_FOUND", 1);
                }
                if( rec.isFlagSet("SPECIES_SAME") ) {
                    getSession().incrementCounter("SPECIES_SAME", 1);
                }
                if( rec.isFlagSet("SPECIES_DIFF") ) {
                    getSession().incrementCounter("SPECIES_DIFF", 1);
                }

                if( rec.isFlagSet("LOAD_UNEXPECTED") ) {
                    getSession().incrementCounter("ERROR_LOAD_UNEXPECTED", 1);
                    dbLog.addLogProp("unknown gene assoc type", "ERROR_LOAD_UNEXPECTED", rec.getRecNo(), PipelineLogger.REC_FLAG);
                }
                else if( rec.isFlagSet("LOAD_SKIP") ) {
                    getSession().incrementCounter("LOAD_SKIP", 1);
                    dbLog.addLogProp("no match against RGD, skipped", "LOAD_SKIP", rec.getRecNo(), PipelineLogger.REC_FLAG);
                }
                else if( rec.isFlagSet("LOAD_MATCH") ) {
                    getSession().incrementCounter("LOAD_MATCH", 1);
                    dbLog.addLogProp("perfect match against RGD, nothing to be done", "LOAD_MATCH", rec.getRecNo(), PipelineLogger.REC_FLAG);
                }
                else if( rec.isFlagSet("LOAD_UPDATE") ) {
                    getSession().incrementCounter("LOAD_UPDATE", 1);
                    dbLog.addLogProp("association subtype has been updated", "LOAD_UPDATE", rec.getRecNo(), PipelineLogger.REC_FLAG);

                    dao.updateAssociation(rec.assoc);
                }
                else if( rec.isFlagSet("LOAD_DIFF_SPECIES") ) {
                    getSession().incrementCounter("LOAD_DIFF_SPECIES", 1);
                }
                else if( rec.isFlagSet("LOAD_ORTHOLOG_SKIP") ) {
                    getSession().incrementCounter("LOAD_ORTHOLOG_SKIP", 1);
                }
                else if( rec.isFlagSet("LOAD_INSERT") ) {
                    getSession().incrementCounter("LOAD_INSERT", 1);
                    dbLog.addLogProp("new assoc inserted", "LOAD_INSERTED", rec.getRecNo(), PipelineLogger.REC_FLAG);

                    dao.insertAssociation(rec.assoc);
                }

                if( !(rec.isFlagSet("LOAD_DIFF_SPECIES") || rec.isFlagSet("LOAD_ORTHOLOG_SKIP")) ) {
                    dbLog.addLogProp("", "", rec.getRecNo(), PipelineLogger.REC_XML, rec.toXml());
                }

                dbLog.writeLogProps(rec.getRecNo());
                dbLog.removeAllLogProps(rec.getRecNo());
            }
        }, "DL", 1, 0);

        manager.getSession().setAllowedExceptions(100);

        // register exception we are *expecting* to happen, so these exceptions
        // won't be reported in the logs
        manager.getSession().registerUserException(new String[]{
                "RGD_ASSOC_UC_IDX",
                "DataIntegrityViolationException",
                "SQLIntegrityConstraintViolationException"});

        manager.run();

        // report associations not handled by the pipeline
        int deletedAssocCount = 0;
        for( Association assoc: assocKeys.values() ) {
            deletedAssocCount += dao.deleteAssociation(assoc);
        }
        manager.getSession().incrementCounter("LOAD_DELETE", deletedAssocCount);


        dbLog.log("Count of records processed", Integer.toString(manager.getSession().getRecordsProcessed(0)), PipelineLog.LOGPROP_RECCOUNT);

        manager.dumpCounters();

        // dump counter statistics
        for( String counter: manager.getSession().getCounters() ) {
            int count = manager.getSession().getCounterValue(counter);
            if( count>0 ) {
                dbLog.log(counter, Integer.toString(count), PipelineLogger.TOTAL);
            }
        }

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

    public class AssocRecord extends PipelineRecord {

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

        String toXml() {
            StringBuilder xml = new StringBuilder();
            xml.append("<?xml version=\"1.0\"?>");
            xml.append("<GeneAssoc>");
            xml.append("<Type>").append(relType).append("</Type>");
            xml.append("<GeneId1>").append(geneId1).append("</GeneId1>");
            xml.append("<RgdId1>");
            appendRgdId(xml, rgdId1);
            xml.append("</RgdId1>");
            xml.append("<GeneId2>").append(geneId2).append("</GeneId2>");
            xml.append("<RgdId2>");
            appendRgdId(xml, rgdId2);
            xml.append("</RgdId2>");
            xml.append("</GeneAssoc>");
            return xml.toString();
        }

        void appendRgdId(StringBuilder xml, int rgdId) {
            if( rgdId==-1 )
                xml.append("GeneId not found in RGD");
            else if( rgdId==-2 )
                xml.append("multiple RGD IDs");
            else if( rgdId==-3 )
                xml.append("withdrawn RGD ID");
            else
                xml.append(rgdId);
        }
    }
}
