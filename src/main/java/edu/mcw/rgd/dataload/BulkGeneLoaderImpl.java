package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogFlagManager;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.XdbManager;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * @author dli
 * Created on Apr 20, 2006
 *
 * heavily rearranged by Marek Tutaj
 */
public class BulkGeneLoaderImpl {

    protected final Logger logger = Logger.getLogger("info");
    protected final Logger logtr = Logger.getLogger("transcripts");
    protected final Logger logtf = Logger.getLogger("transcript_features");

    // settings from AppConfigure.xml
    boolean enableMapPosDeletions; // by default false; if true map positions present in RGD but not found in incoming data will be deleted from db - use with care!

    PipelineLogFlagManager dbFlagManager;

    // map of all gene types found in rgd
    Set<String> rgdGeneTypes = new HashSet<>();

    private CounterPool counters;

    public void update (BulkGene bg) throws Exception {

        updateGeneOnly(bg);

        // insert xdbs
        updateXdbIds(bg.rgdGene.getRgdId(), bg);

        // add new aliases to database
        counters.add("ALIASES_INSERTED", bg.dao.insertAliases(bg.aliases.forInsert));

        if( DataLoadingManager.getInstance().getSpeciesTypeKey()!=8 ) {
            // update fish band map, if the fish map is in RGD, update them, otherwise insert them
            updateMaps(bg);

            // update transcripts and transcript exons
            handleTranscripts(bg);
        }
    }

    public void updateGeneOnly (BulkGene bg) throws Exception {
        Gene rgdGene = bg.rgdGene;

        if (rgdGene==null) {
            logger.error("Rgd Gene not found: GeneId="+bg.getEgId());
            return;
        }

        logger.debug("------Updating genes with rgd id:" + rgdGene.getRgdId());
        if (rgdGene.getRgdId() <= 0) {
            logger.error("The rgd id to be updated is not found");
            return;
        }
        EGDAO dao = bg.dao;

        // update gene type if necessary
        boolean updateGene = false;
        // create nomenclature event
        boolean nomenEvent = false;

        // check if the type is new, if it is new, insert the type first
        String geneTypeLc = checkGeneType(bg.getEgType());
        if( !Utils.stringsAreEqual(rgdGene.getType(), geneTypeLc) ) {
            rgdGene.setType(geneTypeLc);
            updateGene = true;
            logger.debug("gene with rgdid="+rgdGene.getRgdId()+" updated to gene type " + geneTypeLc);
        }

        // handle nomenclature changes
        // update gene symbol and generate nomen event, excluding rat
        String prevSymbol = rgdGene.getSymbol();
        Flags flags = bg.getCustomFlags();
        if (bg.isFlagSet("ORTHO_NOMEN_SYMBOL")) {
            if( bg.getGene().getSpeciesTypeKey()!= SpeciesType.RAT ) {
                rgdGene.setSymbol(bg.getGene().getSymbol());
                updateGene = true;
                nomenEvent = true;
            }

            counters.increment("ORTHO_NOMEN_SYMBOL");
            dbFlagManager.setFlag("ORTHO_NOMEN_SYMBOL", bg.getRecNo());
        }
        if (bg.isFlagSet("ORTHO_NOMEN_SYMBOL_CHANGE_SUPPRESSED_BY_HGNC")) {
            counters.increment("ORTHO_NOMEN_SYMBOL_CHANGE_SUPPRESSED");
            dbFlagManager.setFlag("ORTHO_NOMEN_SYMBOL_CHANGE_SUPPRESSED", bg.getRecNo());
        }

        // update gene name and generate nomen event, excluding rat
        String prevName = rgdGene.getName();
        if (bg.isFlagSet("ORTHO_NOMEN_NAME")) {
            if( bg.getGene().getSpeciesTypeKey()!=SpeciesType.RAT ) {
                rgdGene.setName(bg.getGene().getName());
                updateGene = true;
                nomenEvent = true;
            }
            counters.increment("ORTHO_NOMEN_NAME");
            dbFlagManager.setFlag("ORTHO_NOMEN_NAME", bg.getRecNo());
        }
        if (bg.isFlagSet("ORTHO_NOMEN_NAME_CHANGE_SUPPRESSED_BY_HGNC")) {
            counters.increment("ORTHO_NOMEN_NAME_CHANGE_SUPPRESSED");
            dbFlagManager.setFlag("ORTHO_NOMEN_NAME_CHANGE_SUPPRESSED", bg.getRecNo());
        }

        // update gene description excluding rat
        if (bg.isFlagSet("ORTHO_NOMEN_DESC")) {
            if( bg.getGene().getSpeciesTypeKey()!=SpeciesType.RAT ) {
                rgdGene.setDescription(bg.getGene().getDescription());
                updateGene = true;
            }
            counters.increment("ORTHO_NOMEN_DESC");
            dbFlagManager.setFlag("ORTHO_NOMEN_DESC", bg.getRecNo());
        }

        // update gene refseq status
        if (bg.isFlagSet("GENE_REFSEQ_STATUS")) {
            rgdGene.setRefSeqStatus(bg.getGene().getRefSeqStatus());
            updateGene = true;
            logger.debug("gene with rgdid="+rgdGene.getRgdId()+" updated refseq status to " + rgdGene.getRefSeqStatus());
            counters.increment("GENE_REFSEQ_STATUS_UPDATED");
        }

        // update gene ncbi annot status
        if (flags.getFlagValue().contains(Flags.NCBI_ANNOT_STATUS)) {
            rgdGene.setNcbiAnnotStatus(bg.getGene().getNcbiAnnotStatus());
            updateGene = true;
            logger.debug("gene with rgdid="+rgdGene.getRgdId()+" updated NCBI annot status to " + rgdGene.getNcbiAnnotStatus());
            counters.increment("NCBI_ANNOT_STATUS_UPDATED");
        }

        if( updateGene )
            dao.updateGene(rgdGene);

        if( nomenEvent && bg.getGene().getSpeciesTypeKey()!=SpeciesType.RAT ) {
            // generate nomenclature  event for change of gene symbol/ gene name -- non-rat genes only
            NomenclatureEvent event = new NomenclatureEvent();
            event.setRgdId(rgdGene.getRgdId());
            event.setRefKey("39860");
            event.setEventDate(new Date());
            event.setSymbol(rgdGene.getSymbol());
            event.setName(rgdGene.getName());
            event.setNomenStatusType("APPROVED");
            event.setDesc("Symbol and/or name change");
            event.setOriginalRGDId(rgdGene.getRgdId());
            event.setPreviousSymbol(prevSymbol);
            event.setPreviousName(prevName);
            dao.createNomenEvent(event);
            counters.increment("ORTHO_NOMEN_EVENT");
        }

        //update RGD_IDs table, change last_modified_date
        logger.debug("Updating the last modified date");
        dao.updateLastModifiedDate(rgdGene.getRgdId());
    }

    // insert new gene
    public void load (BulkGene bg) throws Exception {

        Gene newGene = bg.getGene();

        // new gene must have a symbol
        String newGeneSymbol = newGene.getSymbol();
        if (newGeneSymbol== null) {
            logger.error("The gene has no symbol\n");
            return;
        }

        // new gene must have a species given
        int speciesKey = newGene.getSpeciesTypeKey();
        if (speciesKey == SpeciesType.ALL) {
            logger.error("Couldn't get the species for this gene");
            return;
        }

        if( speciesKey!=SpeciesType.RAT ) {
            loadOrthologGene(bg);
            return;
        }

        // check if the gene symbol is already in db
        Gene dupGene = bg.dao.getGeneBySymbolAndSpecies(newGene.getSymbol(), speciesKey);
        if( dupGene!=null ) {
            // if the dupGene is nonactive, change rgd symbol to retired
            // if the dupGene is active, change new gene symbol as NEWGENE_#
            String dupSymbolStatus = bg.dao.getRgdId(dupGene.getRgdId()).getObjectStatus();
            if (dupSymbolStatus.equals("ACTIVE")){
                //if the new gene symbol is duplicate, change gene symbol to NEWGENE_#                
                newGene.setSymbol("NEWGENE_"+dupGene.getRgdId());
            }
            else {
                // change the retired gene symbol to _retired
                updateRetiredGene(dupGene, bg.dao);
            }
        }        
        
        // create RGD ID and load RGD_IDs table first
        logger.debug("------Loading rgd ids table\n");
        String notes = "Created by entrezgene loading";
        RgdId newRgdId = bg.dao.createRgdId(RgdId.OBJECT_KEY_GENES, "ACTIVE", notes, speciesKey);
        int rgdId = newRgdId.getRgdId();
        newGene.setRgdId(rgdId);
        bg.getCustomFlags().setRgdId(rgdId);

        //set gene type, check if the type is new, if it is new, insert the type first
        String geneTypeLc = checkGeneType(bg.getEgType());
        newGene.setType(geneTypeLc);

        // insert genes table
        logger.debug("------Loading genes table, new rgd id is "+ rgdId);
        bg.dao.insertGene(newGene);
        
        // load the nomenEvent table
        NomenclatureEvent event = new NomenclatureEvent();
        event.setRgdId(rgdId);
        event.setRefKey("1724");
        event.setEventDate(new Date());
        event.setSymbol(newGene.getSymbol());
        event.setName(newGene.getName());
        event.setNomenStatusType("PROVISIONAL");
        event.setDesc("Symbol and Name status set to provisional");
        event.setOriginalRGDId(rgdId);
        bg.dao.createNomenEvent(event);
        
        // insert xdbs
        updateXdbIds(rgdId, bg);

        // insert aliases table
        updateAliases(rgdId, bg);

        if( DataLoadingManager.getInstance().getSpeciesTypeKey()!=8 ) {

            // update fish band map, if the fish map is in RGD, update them, otherwise insert them
            updateMaps(rgdId, bg);

            // update transcripts and transcript exons
            handleTranscripts(bg);
        }
    }

    public void updateRetiredGene(Gene gene, EGDAO dao) throws Exception {

        // set gene symbol to _retired
        String prevSymbol = gene.getSymbol();
        gene.setSymbol(prevSymbol+"_retired");
        dao.updateGene(gene);
        logger.debug("One gene symbol make _retired");

        //      update nom events table
        NomenclatureEvent event = new NomenclatureEvent();
        event.setRgdId(gene.getRgdId());
        event.setRefKey("17119");
        event.setEventDate(new Date());
        event.setSymbol(gene.getSymbol());
        event.setPreviousSymbol(prevSymbol);
        event.setName(gene.getName());
        event.setNomenStatusType("APPROVED");
        event.setDesc("Symbol set to symbol_retired");
        event.setOriginalRGDId(gene.getRgdId());
        dao.createNomenEvent(event);
    }

    public void loadOrthologGene(BulkGene bg) throws Exception {

        Gene newGene=bg.getGene();
        int speciesType = newGene.getSpeciesTypeKey();

        // check gene type table, if gene type doesn't exist, insert the new gene type
        String geneTypeLc = checkGeneType(bg.getEgType());
        newGene.setType(geneTypeLc);

        // create RGD ID and load RGD_IDs table first
        logger.debug("------Loading rgd ids table\n");
        RgdId newRgdId = bg.dao.createRgdId(RgdId.OBJECT_KEY_GENES, "ACTIVE", "Created by entrezgene loading", speciesType);
        int rgdId = newRgdId.getRgdId();
        bg.getCustomFlags().setRgdId(rgdId);

        // insert genes table
        logger.debug("------Loading genes table, new rgd id is "+ rgdId);
        newGene.setRgdId(rgdId);       
        bg.dao.insertGene(newGene);

        // insert xdbs
        updateXdbIds(rgdId, bg);

        // insert aliases table
        updateAliases(rgdId, bg);

        if( DataLoadingManager.getInstance().getSpeciesTypeKey()!=8 ) {

            // update maps
            updateMaps(rgdId, bg);

            // update transcripts and transcript exons
            handleTranscripts(bg);
        }
    }



    void updateMaps(BulkGene bg) throws Exception {
        logger.debug("Updating map data");
        bg.genePositions.syncMapData(bg, null, getDbFlagManager(), "GENE", true, counters);
    }

    // for orthologs
    void updateMaps(int rgdId, BulkGene bg) throws Exception {

        bg.genePositions.updateRgdId(rgdId);
        updateMaps(bg);
    }

    // update gene type
    // check if the type is new, if it is new, insert the type first
    String checkGeneType(String egGeneType) throws Exception {
        EGDAO egDAO = EGDAO.getInstance();
        if (egGeneType!=null) {
            String geneTypeLc = egGeneType.toLowerCase();

            // check if the gene type in on our map of known rgd genetypes
            boolean typeInRgd = rgdGeneTypes.contains(geneTypeLc);
            if( !typeInRgd ) {
                // gene type is not in our map: get it from database
                typeInRgd = egDAO.existsGeneType(geneTypeLc);
                if( typeInRgd )
                    // gene type found in rgd database: update our gene type map
                    rgdGeneTypes.add(geneTypeLc);
            }

            if (!typeInRgd) {
                // the type is not in RGD
                egDAO.createGeneType(geneTypeLc, geneTypeLc, geneTypeLc);
            }
            return geneTypeLc;
        }
        return null;
    }

    void updateXdbIds(int rgdId, BulkGene bg) throws Exception {

        // are any xdb ids to be inserted?
        List<XdbId> xdbs = bg.getXdbIds();
        if (xdbs!=null ) {

            GeneXdbIds geneXdbIds = new GeneXdbIds();

            // ensure all xdb ids have valid rgd-id
            for( XdbId id: xdbs ) {
                id.setRgdId(rgdId);
                id.setSrcPipeline(XdbManager.EG_PIPELINE);

                geneXdbIds.addIncoming(id);
            }

            geneXdbIds.qc(rgdId, counters);
        }
    }

    void updateAliases(int rgdId, BulkGene bg) throws Exception {

        // update rgd id for the aliases to add
        for (Alias alias: bg.aliases.forInsert) {
            alias.setRgdId(rgdId);
        }

        // add new aliases to database
        counters.add("ALIASES_INSERTED", bg.dao.insertAliases(bg.aliases.forInsert));
    }


    void handleTranscripts(BulkGene bg) throws Exception {

        Flags flag = bg.getCustomFlags();
        int rgdid = flag.getRgdId();
        logger.debug("transcript checking: EG_ID="+bg.getEgId()+", GENE_RGD_ID="+rgdid);

        // see if all incoming transcripts are in rgd; if not add them
        updateTranscripts(bg);

        updateTranscriptPositions(bg);
        updateTranscriptFeatures(bg);
        updateTranscriptXdbIds(bg);

        logger.debug("========= TRANSCRIPTS DONE: EG_ID="+bg.getEgId()+", GENE_RGD_ID="+rgdid);
    }

    void updateTranscripts(BulkGene bg) throws Exception {

        int rgdId = bg.getCustomFlags().getRgdId();

        List<Transcript> obsoleteInRgdTranscripts = new ArrayList<>(bg.rgdTranscripts);

        // see if all incoming transcripts are in rgd
        for( TranscriptInfo ti: bg.transcripts ) {
            ti.setGeneRgdId(rgdId);

            boolean transcriptInRgd = false;
            for( Transcript tr: bg.rgdTranscripts ) {
                if( Utils.stringsAreEqual(tr.getAccId(), ti.getAccId()) ) {
                    obsoleteInRgdTranscripts.remove(tr);

                    // set transcript rgd id for incoming transcript
                    ti.setRgdId(tr.getRgdId());

                    // update transcript properties if needed
                    if( !tr.equals(ti) ) {
                        // properties differ between incoming transcript and rgd
                        // log the changed properties
                        logtr.info("OLD_TRANSCRIPT_PROPERTIES: "+tr.dump("|"));

                        if( !Utils.stringsAreEqual(tr.getRefSeqStatus(), ti.getRefSeqStatus()) ) {
                            getDbFlagManager().setFlag("TRANSCRIPT_REFSEQ_STATUS_CHANGED", bg.getRecNo());
                            counters.increment("TRANSCRIPT_REFSEQ_STATUS_CHANGED");
                            tr.setRefSeqStatus(ti.getRefSeqStatus());
                        }

                        if( tr.isNonCoding()!=ti.isNonCoding() ) {
                            getDbFlagManager().setFlag("TRANSCRIPT_CODING_STATUS_CHANGED", bg.getRecNo());
                            counters.increment("TRANSCRIPT_CODING_STATUS_CHANGED");
                            tr.setNonCoding(ti.isNonCoding());
                        }

                        if( !Utils.stringsAreEqual(tr.getPeptideLabel(), ti.getPeptideLabel()) ) {
                            getDbFlagManager().setFlag("TRANSCRIPT_PEPTIDE_LABEL_CHANGED", bg.getRecNo());
                            counters.increment("TRANSCRIPT_PEPTIDE_LABEL_CHANGED");
                            tr.setPeptideLabel(ti.getPeptideLabel());
                        }

                        if( !Utils.stringsAreEqual(tr.getProteinAccId(), ti.getProteinAccId()) ) {
                            getDbFlagManager().setFlag("TRANSCRIPT_PROTEIN_ACC_ID_CHANGED", bg.getRecNo());
                            counters.increment("TRANSCRIPT_PROTEIN_ACC_ID_CHANGED");
                            tr.setProteinAccId(ti.getProteinAccId());
                        }

                        if( bg.dao.updateTranscript(tr)!=0 ) {
                            counters.increment("TRANSCRIPT_UPDATED");
                        }

                        logtr.info("NEW_TRANSCRIPT_PROPERTIES: "+tr.dump("|"));
                    }

                    transcriptInRgd = true;
                    break;
                }
            }

            if( !transcriptInRgd ) {
                Transcript newTr = new Transcript();
                newTr.setGeneRgdId(ti.getGeneRgdId());
                newTr.setNonCoding(ti.isNonCoding());
                newTr.setAccId(ti.getAccId());
                newTr.setRefSeqStatus(ti.getRefSeqStatus());
                newTr.setProteinAccId(ti.getProteinAccId());
                newTr.setPeptideLabel(ti.getPeptideLabel());

                bg.dao.createTranscript(newTr, bg.getGene().getSpeciesTypeKey());
                ti.setRgdId(newTr.getRgdId());
                bg.rgdTranscripts.add(newTr); // we new transcript now in rgd!
                logtr.info("TRANSCRIPT_INSERTED: "+newTr.dump("|"));

                // increment counter of transcripts added
                counters.increment("TRANSCRIPTS_INSERTED");
            }
            else {
                counters.increment("TRANSCRIPTS_MATCHED");
            }
        }

        // detach in-rgd transcripts that are no longer in the incoming list
        if( !obsoleteInRgdTranscripts.isEmpty() ) {

            // but only if there are incoming transcripts available
            // we want to avoid detaching all transcript for given gene when
            // f.e. the gene is no longer on the current assembly
            if( bg.transcripts.isEmpty() ) {
                counters.increment("TRANSCRIPT_DETACH_FROM_GENE_SUPPRESSED");
                getDbFlagManager().setFlag("TRANSCRIPT_DETACH_FROM_GENE_SUPPRESSED", bg.getRecNo());
            } else {
                for( Transcript tr: obsoleteInRgdTranscripts ) {
                    bg.dao.detachTranscriptFromGene(tr);
                    logtr.info("TRANSCRIPT_DETACHED_FROM_GENE: "+tr.dump("|"));

                    // increment counter of transcripts detached
                    counters.increment("TRANSCRIPTS_DETACHED");
                }
                getDbFlagManager().setFlag("TRANSCRIPT_DETACHED_FROM_GENE", bg.getRecNo());
            }
        }
    }

    void updateTranscriptPositions(BulkGene bg) throws Exception {

        final Logger logger = Logger.getLogger("transcript_positions");

        GenePositions positions = new GenePositions(bg.dao);

        // populate table with positions for incoming transcripts
        for( TranscriptInfo ti: bg.transcripts ) {
            for( TranscriptLocus locus: ti.getLoci() ) {
                TranscriptFeature tf = locus.getTranscriptCoords();
                tf.setRgdId(ti.getRgdId());
                positions.addMapData(tf);
            }
        }

        // populate table with positions for RGD transcripts
        List<MapData> trPosInRgd = new ArrayList<MapData>();
        for( Transcript tr: bg.rgdTranscripts ) {
            trPosInRgd.addAll(tr.getGenomicPositions());
        }
        positions.setRgdMapData(trPosInRgd, null);


        // now synchronize transcript positions between incoming data and RGD
        positions.qcMapData(bg, logger);
        positions.syncMapData(bg, logger, getDbFlagManager(), "TRANSCRIPT", true, counters);

        positions.deleteOverlappingPositionsMarkedForDelete(bg, getDbFlagManager(), counters);
    }

    void updateTranscriptFeatures(BulkGene bg) throws Exception {

        List<TranscriptFeature> matchingRgdFeatures = new ArrayList<TranscriptFeature>(); // rgd features that do match the incoming features
        List<TranscriptFeature> rgdFeatures; // features in rgd to be processed
        // only those that should be sync-ed -- if a feature is in RGD but its map key refers to old assembly,
        // it won't be added to this array to avoid being deleted from database

        // determine list of map keys to be processed
        rgdFeatures = getFeaturesMatchingLocus(bg.rgdFeatures, bg.transcripts);
        List<TranscriptFeature> inRgdFeatures = new ArrayList<>(rgdFeatures);

        // now determine matching features
        for( TranscriptInfo tr: bg.transcripts ) {

            for( TranscriptLocus locus: tr.getLoci() ) {

                // try to match all genomic features, like exons
                List<TranscriptFeature> features = new LinkedList<>(locus.getCoords());
                if( locus.getUtr3()!=null ) {
                    features.add(locus.getUtr3());
                }
                if( locus.getUtr5()!=null ) {
                    features.add(locus.getUtr5());
                }

                List<TranscriptFeature> rgdFeaturesInTranscriptLocus = getRgdFeaturesInTranscriptLocus(inRgdFeatures, locus, bg, tr.getRgdId());

                for( TranscriptFeature tf: features ) {
                    tf.setTranscriptRgdId(tr.getRgdId());

                    // match the feature with rgd feature
                    TranscriptFeature matchedFeature = matchFeature(tf, rgdFeatures, bg);
                    if( matchedFeature==null ) {
                        continue;
                    }
                    matchingRgdFeatures.add(matchedFeature);
                    if( matchedFeature.getFeatureType()== TranscriptFeature.FeatureType.EXON ) {
                        counters.increment("EXONS_MATCHED");
                    } else if( matchedFeature.getFeatureType()== TranscriptFeature.FeatureType.UTR3 ) {
                        counters.increment("UTRS_MATCHED");
                    } else if( matchedFeature.getFeatureType()== TranscriptFeature.FeatureType.UTR5 ) {
                        counters.increment("UTRS_MATCHED");
                    }

                    // remove the matching feature from in-rgd features
                    Iterator<TranscriptFeature> it = rgdFeaturesInTranscriptLocus.iterator();
                    while( it.hasNext() ) {
                        TranscriptFeature tfInRgd = it.next();
                        if( tfInRgd.getTranscriptRgdId()==matchedFeature.getTranscriptRgdId() &&
                                tfInRgd.getFeatureType()==matchedFeature.getFeatureType() &&
                                tfInRgd.equalsByGenomicCoords(matchedFeature) ) {
                            // match
                            it.remove();
                        }
                    }
                }

                // any stale in-rgd features in this locus must be unlinked
                unlinkFeatures(rgdFeaturesInTranscriptLocus, bg);
            }
        }

        // determine unmatching features in rgd to be unlinked
        rgdFeatures.removeAll(matchingRgdFeatures);
        unlinkFeatures(rgdFeatures, bg);
    }

    void unlinkFeatures(List<TranscriptFeature> rgdFeatures, BulkGene bg) throws Exception {

        if( rgdFeatures.isEmpty() ) {
            return;
        }

        List<TranscriptFeature> rgdFeaturesToUnlink = new ArrayList<TranscriptFeature>();
        for( TranscriptFeature f: rgdFeatures ) {
            logtf.info("Feature to unlink: GeneId="+bg.getEgId()+" RgdId="+bg.gene.getRgdId()+" Gene="+bg.gene.getSymbol()+", "+f);
        }
        rgdFeaturesToUnlink.addAll(rgdFeatures);

        // determine which transcript features have not been matched and remove them from database
        int rowsAffected = 0;
        for( TranscriptFeature tf: rgdFeaturesToUnlink ) {
            rowsAffected += bg.dao.unlinkFeature(tf.getRgdId());

            counters.increment(
                    tf.getFeatureType()== TranscriptFeature.FeatureType.EXON
                            ? "EXONS_UNLINKED"
                            : "UTRS_UNLINKED");
        }
        if( rowsAffected>0 ) {
            logger.debug("===== UNLINK FEATURES: rows affected " + rowsAffected);
        }
    }

    List<TranscriptFeature> getRgdFeaturesInTranscriptLocus(List<TranscriptFeature> rgdFeatures, TranscriptLocus locus, BulkGene bg, int trRgdId) {

        // determine gene position matching transcript locus
        MapData md = bg.genePositions.getPositionMatchingLocus(locus.getTranscriptCoords());
        List<TranscriptFeature> result = new ArrayList<>();
        for( TranscriptFeature tf: rgdFeatures ) {
            if( tf.getTranscriptRgdId()==trRgdId && GenePositions.positionsOverlap(tf, md) ) {
                result.add(tf);
            }
        }
        return result;
    }

    /**
     * from rgdFeatures, get a subset matching maps from incoming data
     * @param rgdFeatures rgdFeatures
     * @param trs list of incoming transcripts -- map keys are examined
     * @return list of rgdFeatures matching a list of allowable map keys
     */
    List<TranscriptFeature> getFeaturesMatchingLocus(List<TranscriptFeature> rgdFeatures, List<TranscriptInfo> trs) {

        // determine list of map keys to be processed
        Set<Integer> mapKeys = new HashSet<Integer>();
        for( TranscriptInfo tr: trs ) {
            for( TranscriptLocus locus: tr.getLoci() ) {
                mapKeys.add(locus.getTranscriptCoords().getMapKey());
            }
        }

        List<TranscriptFeature> results = new ArrayList<TranscriptFeature>();
        for( TranscriptFeature tf: rgdFeatures ) {
            if( mapKeys.contains(tf.getMapKey()) ) {
                results.add(tf);
            }
        }
        return results;
    }

    TranscriptFeature matchFeature(TranscriptFeature tfIncoming, List<TranscriptFeature> rgdTFList, BulkGene bg) throws Exception {
        // lookup for features already attached to this transcript in RGD
        for( TranscriptFeature tfInRgd: rgdTFList ) {
            if( tfIncoming.equalsByGenomicCoords(tfInRgd) && tfIncoming.getTranscriptRgdId()==tfInRgd.getTranscriptRgdId() ) {
                // incoming feature matches the rgd feature: same transcript, map, chromosome, strand, start and stop position
                tfIncoming.setRgdId(tfInRgd.getRgdId()); // reuse feature_rgd_id
                
                return tfIncoming;
            }
        }

        // lookup for features to be shared (bound to another transcript already)
        for( TranscriptFeature tfInRgd: rgdTFList ) {
            if( tfIncoming.equalsByGenomicCoords(tfInRgd) ) {
                // incoming feature matches the rgd feature: same map, chromosome, strand, start and stop position
                tfIncoming.setRgdId(tfInRgd.getRgdId()); // reuse feature_rgd_id
                // create a binding between the feature and the current transcript
                bg.dao.bindFeatureToTranscript(tfIncoming.getTranscriptRgdId(), tfInRgd.getRgdId());
                return tfIncoming;
            }
        }

        // unmatching feature -- create a new transcript feature object
        bg.dao.createFeature(tfIncoming, bg.getGene().getSpeciesTypeKey());
        rgdTFList.add(tfIncoming); // we new transcript now in rgd!
        counters.increment(
                tfIncoming.getFeatureType()== TranscriptFeature.FeatureType.EXON
                    ? "EXONS_INSERTED"
                    : "UTRS_INSERTED");
        logtf.info("FEATURE INSERT "+tfIncoming.toString());
        return tfIncoming;
    }

    void updateTranscriptXdbIds(BulkGene bg) throws Exception {

        for( Transcript tr: bg.rgdTranscripts ) {
            bg.transcriptXdbIds.qc(tr, bg);
            bg.transcriptXdbIds.load(bg, counters);
        }
    }

    public boolean isEnableMapPosDeletions() {
        return enableMapPosDeletions;
    }

    public void setEnableMapPosDeletions(boolean enableMapPosDeletions) {
        this.enableMapPosDeletions = enableMapPosDeletions;
    }

    public PipelineLogFlagManager getDbFlagManager() {
        return dbFlagManager;
    }

    public void setDbFlagManager(PipelineLogFlagManager dbFlagManager) {
        this.dbFlagManager = dbFlagManager;
    }

    public CounterPool getCounters() {
        return counters;
    }

    public void setCounters(CounterPool counters) {
        this.counters = counters;
    }
}
