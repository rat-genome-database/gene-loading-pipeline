package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.*;
import org.apache.log4j.Logger;

import java.sql.Timestamp;
import java.util.*;
import java.util.Map;

/**
 * <ol>
 * <li> Check if the gene is in RGD by matching Entrezgene ID
 * <li> Check if the new entrezgene has rgd id
 * <li> Check if the new entrezgene rgd id match rgd gene rgd id
 * </ol>
 */
public class QualityCheckBulkGene  {
    XdbManager xdbManager = new XdbManager();
    TranscriptDAO transcriptDAO = new TranscriptDAO();
    PipelineLogger dbLog = PipelineLogger.getInstance();
    PipelineLogFlagManager dbFlagManager;

    // autoloaded from AppConfig.xml
    List egType;
    List excludeEgType;
    private String version;

    // genomic assemblies to be handled by the pipeline, specific to the species
    // (key: assembly name; value: map key)
    Map<String, String> genomicAssemblies;
    // map keys to be processed by the pipelines
    // (derived from the list of genomic assemblies)
    // NOTE: this avoids deletions of other older assembly map positions for genes and transcripts
    Set<Integer> validMapKeys;

    protected final Logger logger = Logger.getLogger("process");
    protected final Logger logInactive = Logger.getLogger("inactive");
    protected final Logger logNomen = Logger.getLogger("nomen"); // tracks changes in gene names, symbols and descriptions
    protected final Logger logSymbol = Logger.getLogger("symbol"); // tracks changes in gene symbols, user friendly

    // do pipeline checks and see
    public Flags check (BulkGene fileGene) throws Exception {

        // check if the new entrezGene has an entrezGene ID; do logs
        Flags flag = checkEgId(fileGene);
        if( flag!=null )
            return flag;

        int speciesTypeKey = fileGene.getGene().getSpeciesTypeKey();
        if( speciesTypeKey==SpeciesType.ALL ) {
            // the species had not been recognized during xml parsing;
            // check if we have species-type-key override
            // if yes, use the species type key override
            int speciesTypeKeyOverride = DataLoadingManager.getInstance().getSpeciesTypeKey();
            if( speciesTypeKeyOverride!=SpeciesType.ALL ) {
                fileGene.getGene().setSpeciesTypeKey(speciesTypeKeyOverride);
                speciesTypeKey = fileGene.getGene().getSpeciesTypeKey();
                logger.debug("###########incoming EG species (override): " + SpeciesType.getTaxonomicName(speciesTypeKey));
            }
        }

        switch (speciesTypeKey) {
            case SpeciesType.RAT: {
                return ratCheck(fileGene);
            }
            case SpeciesType.HUMAN: {
                XdbId xdb = new XdbId();
                xdb.setXdbKey(XdbId.XDB_KEY_HGNC);
                return orthologsCheck(fileGene, xdb);
            }
            case SpeciesType.MOUSE: {
                XdbId xdb = new XdbId();
                xdb.setXdbKey(XdbId.XDB_KEY_MGD);
                return orthologsCheck(fileGene, xdb);
            }
            default: {
                XdbId xdb = new XdbId();
                return orthologsCheck(fileGene, xdb);
            }
        }
    }
    
    public Flags ratCheck (BulkGene bulkGene) throws Exception {

        Gene ngene = bulkGene.getGene(); // gene initialized from EntrezGene record
        Flags flag = new Flags();
        EGDAO dao = EGDAO.getInstance();

        // check if the eg type is one of the types we want
        checkAllowedTypes(bulkGene);

        //check if the egID match any egID in RGD
        //--------------------need to check status----------------------------
        List<Gene> genesByEgId = dao.getGenesByEGID(bulkGene.getEgId());
        List<Gene> activeRgdGene = new ArrayList<>();
        
        int geneTypeC=0;
        String info="";
        // get matched active gene and count the active genes of type "gene"
        for( Gene rgdGeneTmp: genesByEgId ) {
            // get the rgd gene status, and put the active gene into a list
            String rgdGeneStatus = dao.getRgdId(rgdGeneTmp.getRgdId()).getObjectStatus();
            if (rgdGeneStatus.equals("ACTIVE")) {
                activeRgdGene.add(rgdGeneTmp);
                String rgdGeneType=rgdGeneTmp.getType();
                if (rgdGeneType!=null ) {
                    bulkGene.rgdGene = rgdGeneTmp;
                    geneTypeC ++;
                }
                info += rgdGeneTmp.getRgdId()+",";
            }
        }
        
        if (activeRgdGene.size()>1) {
            if (geneTypeC == 1) {
                //one EG ID has multiple genes in RGD, only one is type "gene"
                // proceed to load the gene
                logger.debug("incoming EG ID has multiple genes in RGD, only one is type (gene)");
                flag = new Flags(Flags.EGINRGD_VAR);
                flag.setRelatedInfo("Multiple genes in RGD: "+ info);
                dbLog.addLogProp(flag.getRelatedInfo(), "EG_IN_RGD_VAR", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            }
            else if (geneTypeC > 1 || geneTypeC==0) {
                // one EG ID has multiple genes in RGD, more than one gene has type "gene" or no gene has type "gene"
                // write to error log file
                logger.debug("incoming EG is associated multiple genes in RGD");
                flag = new Flags(Flags.MULTIGENES, Flags.ERROR);
                flag.setRelatedInfo("Associated with multiple RGD genes: "+ info);
                dbLog.addLogProp(flag.getRelatedInfo(), "MULTIGENES", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                return flag;
            }
        }  
        else if (activeRgdGene.size()== 1){
            logger.debug("incoming EG is associated with one gene in RGD");
            flag = new Flags(Flags.EGINRGD);
            dbLog.addLogProp("incoming EG is associated with one gene in RGD", "EG_IN_RGD", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            bulkGene.rgdGene = activeRgdGene.get(0);
        }
        else {
            // the egid doesn't match active egid in rgd            
           
            if (ngene.getRgdId()>0) {
                // The eg id is not in RGD, the new eg record have a rgdid
                // update eg id in RGD, put the EGID in a holding array
                logger.debug("eg id is not in RGD, incoming eg record has a rgdid");
                try {
                    bulkGene.rgdGene = dao.getGene(ngene.getRgdId());
                }
                catch( GeneDAO.GeneDAOException e ) {
                    // there is no such RGD_ID in RGD database!
                    logger.debug("eg id is not in RGD, incoming eg record has an rgd-id, but there is no such RGD_ID in rgd database!");
                    // check if the gene symbol is in RGD
                    flag = new Flags(Flags.NEWGENE, Flags.INSERT);
                    dbLog.addLogProp("eg id is not in RGD, incoming eg record has an rgd-id, but there is no such RGD_ID in rgd database!", "NEW_GENE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                    return flag;
                }

                String rgdGeneStatus = dao.getRgdId(ngene.getRgdId()).getObjectStatus();
                // check if the rgd gene is active
                if (!rgdGeneStatus.equals("ACTIVE")) {
                    // if the rgd gene is non active gene
                    // the rgd gene either has different eg id or has no eg id
                    dbLog.addLogProp("eg id is not in RGD, incoming eg record has a "+rgdGeneStatus+" rgdid", "NON_ACTIVE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                    logInactive.info("RGDID="+ngene.getRgdId()+", EGID="+bulkGene.getEgId());
                    return new Flags(Flags.NONACTIVE, Flags.ERROR);
                }
                else {
                    // the rgd gene is active, has different eg id or has no eg id
                    flag = new Flags(Flags.DIFFEGID);
                    //flag.setLoadStatus(Flags.ERROR);
                    flag.setRelatedInfo("incoming EG has RGD ID, EG ID not in RGD");
                    dbLog.addLogProp(flag.getRelatedInfo(), "DIFF_EG_ID", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                }
            }
            else {
                if (genesByEgId.size()>0) {
                    // the eg id is not active in RGD and new eg record has no rgdid
                    // write to error log
                    flag = new Flags(Flags.NONACTIVE, Flags.ERROR);
                    dbLog.addLogProp("eg id is not active in RGD and incoming eg record has no rgdid", "NON_ACTIVE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                    for( RgdId rgdId: dao.getRGDIdsByXdbId(XdbId.XDB_KEY_ENTREZGENE, bulkGene.getEgId()) ) {
                        logInactive.info("RGDID="+rgdId.getRgdId()+", EGID="+bulkGene.getEgId());
                    }
                    return flag;

                } else {
                    //the eg id is not in RGD and the new eg record doesn't have a rgdid
                    logger.debug("eg id is not in RGD, incoming eg record doesn't have a rgdid");
                    // check if the gene symbol is in RGD
                    flag = new Flags(Flags.NEWGENE, Flags.INSERT);
                    dbLog.addLogProp("eg id is not in RGD, incoming eg record doesn't have a rgdid", "NEW_GENE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                    return flag;
                }
            }
        }
        //to here, entrezgene ID match active gene or 
        //the entrezgene ID dosn't match any eg id, but new record has rgd id and that rgd id is active
        // the matched rgd gene is definitely active now
        // the flag value is either EGINRGD or EGINRGD_VAR
        //check if the rgd id match        
        flag.setRgdId(bulkGene.rgdGene.getRgdId());
        flag.setRelatedInfo("");
      
        if( ngene.getRgdId()!=0 && bulkGene.rgdGene.getRgdId()!=ngene.getRgdId()) {
            // incoming rgd id doesn't match the rgd gene -- that could be because
            // incoming rgd id was replaced in RGD by other gene
            int newActiveGeneRgdId = dao.getReplacedGeneByRgdId(ngene.getRgdId());
            if( newActiveGeneRgdId!=bulkGene.rgdGene.getRgdId() ) {
                logger.debug("EGID match, rgd id doesn't match");
                flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_DIFFRGDID);
                flag.setLoadStatus(Flags.ERROR);
                dbLog.addLogProp("eg id match, rgd id doesn't match", "EG_IN_RGD_DIFF_RGD_ID", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                return flag;
            }
        }

        // rgd id match or new record does not have rgd id
        // check if the rgd record has sequence
        List<XdbId> rgdSeqXdbs = dao.getXdbIdsByRgdId(XdbId.XDB_KEY_GENEBANKNU, bulkGene.rgdGene.getRgdId());
        rgdSeqXdbs = xdbManager.removeGeneBankNucleotides(rgdSeqXdbs);
        List newSeqXdbs = bulkGene.getXdbIdsByXdbKey(XdbId.XDB_KEY_GENEBANKNU);
        //logger.debug("the rgd xdb size: "+ rgdSeqXdbs.size());
        if (rgdSeqXdbs.size()==0 && newSeqXdbs.size()>0) {
            // the eg id and rgd id match, rgd has no sequence, new eg record has sequence
            // remove the EG ID from rgd gene
            // make new RGD record for the EG record
            // write to log file mapped_gene_upgrade
            logger.debug("eg id and rgd id match, the rgd record has no sequence");
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp("eg id and rgd id match, the rgd record has no sequence", "UPGRADE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            // note: UPGRADE step obsolete -- there should be no logic with _mapped
            return flag;
        }
        else if (rgdSeqXdbs.size()==0 && newSeqXdbs.size()==0) {
            logger.debug("eg id and rgd id match, neither rgd record nor incoming eg record has sequence");
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp("eg id and rgd id match, neither rgd record nor incoming eg record has sequence", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }
        else {
            logger.debug("eg id and rgd id match, rgd record has sequence");
            //ArrayList<RgdAccXdb> newSeqXdbs=fileGene.getRgdAccXdbByName(RgdAccXdb.GBNU);

            if (oneXdbMatch(rgdSeqXdbs, newSeqXdbs)) {
                //the eg id and rgd id match, at least one sequence match
                // proceed to update
                flag.setLoadStatus(Flags.UPDATE);
                // check the new eg and rgd gene type
                if (bulkGene.getEgType().contains("pseudotype")) {
                    if (bulkGene.rgdGene.getType().contains("pseudotype")) {
                        logger.debug("eg id and rgd id match, at least one sequence match, the type of both incoming eg and rgd gene is pseudo");
                        // proceed to update
                        dbLog.addLogProp("eg id and rgd id match, at least one sequence match, the type of both incoming eg and rgd gene is pseudo", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                    }
                    else {

                        if (bulkGene.rgdGene.getType().contains("predicted")) {
                            logger.debug("eg id and rgd id match, at least one sequence match, eg type is pseudo, rgd type is predicted");
                            // proceed to update
                            dbLog.addLogProp("eg id and rgd id match, at least one sequence match, eg type is pseudo, rgd type is predicted", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);

                        }
                        else {
                            // change gene type, do update
                            flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_DIFFTYPE);
                            dbLog.addLogProp("eg id and rgd id match, at least one sequence match, eg type is pseudo, rgd type is neither pseudo nor predicted", "EG_IN_RGD_DIFF_TYPE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                        }
                    }
                }
                else {
                    logger.debug("eg id and rgd id match, at least one sequence match, eg type is not pseudo");
                    dbLog.addLogProp("eg id and rgd id match, at least one sequence match, eg type is pseudo, incoming eg type is not pseudo", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                }

                return flag;
            }
            else {
                // there is not a single match between sequences in incoming data vs RGD
                //
                // is there a single XM or XR or NG sequence in RGD?
                boolean isXmXrNg = hasSequence(rgdSeqXdbs, new String[]{"XM_","XR_","NG_"});

                // there is XM or XR or NG sequence in RGD
                if (isXmXrNg) {
                    return handleXmXrNgSeq(flag, bulkGene);
                }
                else {
                    // the eg id and rgd id match, sequence doesn't match, no sequence is XM or XR or NG
                    //
                    handleMissingSequence(flag, bulkGene, rgdSeqXdbs, newSeqXdbs);
                    return flag;
                }
            }
        }
    }

    // return true if there is at least one of the specified sequences given in 'seqs' array on the list
    boolean hasSequence(List<XdbId> xdbs, String[] seqs) {

        for (XdbId rgdXdb: xdbs) {
            if (rgdXdb.getAccId()!=null) {
                for( String seq: seqs ) {
                    if (rgdXdb.getAccId().startsWith(seq) ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    // return true if all of the specified sequences given in 'seqs' array are on the list
    boolean hasAllSequences(List<XdbId> xdbs, String[] seqs) {

        for (XdbId rgdXdb: xdbs) {
            if (rgdXdb.getAccId()!=null) {
                boolean isMatch = false;
                for( String seq: seqs ) {
                    if (rgdXdb.getAccId().startsWith(seq) ) {
                        isMatch = true;
                        break;
                    }
                }
                if( !isMatch )
                    return false;
            }
        }
        return true;
    }

    private void handleMissingSequence(Flags flag, BulkGene bulkGene, List<XdbId> rgdSeqXdbs, List<XdbId> newSeqXdbs) throws Exception {

        // the eg id and rgd id match, sequence doesn't match, no sequence in RGD is XM or XR or NG
        String changedSeqIds = "RGD SEQs:"+rgdAccXdbToString(rgdSeqXdbs)+",NEW SEQs:"+rgdAccXdbToString(newSeqXdbs);
        logger.debug("eg id and rgd id match, sequence doesn't match, no sequence in RGD is XM or XR or NG");

        // IF the gene is pseudo AND there is no incoming sequence present, but THERE are sequences in RGD,
        // THEN flag it as PSEUDO_GENE_WITH_NO_SEQ and proceed to gene UPDATE
        if (bulkGene.getEgType().contains("pseudo") && !rgdSeqXdbs.isEmpty() && newSeqXdbs.isEmpty()) {

            // proceed to update
            flag.setLoadStatus(Flags.UPDATE);

            String msg = "eg id and rgd id match, pseudo gene, sequence in RGD, no incoming sequence";
            logger.debug(msg);
            flag.setFlagValue(flag.getFlagValue()+","+ Flags.PSEUDO_WITH_NOSEQ);
            dbLog.addLogProp(msg+"; "+changedSeqIds, Flags.PSEUDO_WITH_NOSEQ, bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
        }
        else if( hasSequence(newSeqXdbs, new String[]{"XM_","XR_","NG_","NM_","NR_"}) ) {
            // EG ID and rgd id match, sequence does not match, no sequence in RGD is XM or XR or NG
            // but there is XM or XR or NG or NM or NR sequence in incoming data

            flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_NEWSEQS);
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp(changedSeqIds+"<br>eg id and rgd id match, sequence doesn't match, no sequence in RGD is XM or XR or NG, incoming sequence is XM or XR or NG or NR or NM", "EG_IN_RGD_NEW_SEQS", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            // print the accession ID
            flag.setRelatedInfo(flag.getRelatedInfo()+", "+changedSeqIds);
        }
        else
        // EG ID and rgd id match, sequence does not match, no sequence in RGD is XM or XR or NG or NM or NR
        // and there is no XM or XR or NG sequence in incoming data
        // but both the incoming data and matching gene RGD do contain AA sequences
        if( hasSequence(newSeqXdbs, new String[]{"AABR", "AAHX"}) && hasAllSequences(rgdSeqXdbs, new String[]{"AABR", "AAHX"}) ) {

            flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_AASEQS);
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp(changedSeqIds+"<br>eg id and rgd id match, sequence doesn't match, new AA sequences", "EG_IN_RGD_AA_SEQS", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            // print the accession ID
            flag.setRelatedInfo(flag.getRelatedInfo()+", "+changedSeqIds);

        } else {
            flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_DIFFSEQ);
            flag.setLoadStatus(Flags.ERROR);
            dbLog.addLogProp(changedSeqIds+"<br>eg id and rgd id match, sequence doesn't match, no sequence is XM or XR or NG or NM or NR", "EG_IN_RGD_DIFF_SEQ", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            // print the accession ID
            flag.setRelatedInfo(flag.getRelatedInfo()+", "+changedSeqIds);
        }
    }

    Flags handleXmXrNgSeq(Flags flag, BulkGene bulkGene) {

        // the eg id and rgd id match, sequence doesn't match, has XM or XR or NG sequence
        logger.debug("eg id and rgd id match, sequence doesn't match, has XM or XR or NG sequence");
        flag.setLoadStatus(Flags.UPDATE);
        dbLog.addLogProp("eg id and rgd id match, sequence doesn't match, has XM or XR or NG sequence", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);

        if (bulkGene.getEgType().contains("pseudotype")) {
            if (bulkGene.rgdGene.getType().contains("pseudotype")) {
                logger.debug("eg id and rgd id match, sequence doesn't match, has xm sequence, the type of both incoming eg and rgd gene is pseudo");
                // proceed to update
                dbLog.addLogProp("... type of both incoming eg and rgd gene is pseudo", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            }
            else {

                if (bulkGene.rgdGene.getType().contains("predicted")) {
                    logger.debug("eg id and rgd id match, sequence doesn't match, has xm sequence, eg type is pseudo, rgd type is predicted");
                    // proceed to update
                    dbLog.addLogProp("... eg type is pseudo, rgd type is predicted", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                }
                else {
                    // change gene type, do update
                    flag.setFlagValue(flag.getFlagValue()+","+ Flags.EGINRGD_DIFFTYPE);
                    flag.setRelatedInfo(flag.getRelatedInfo()+",Changed gene type in RGD");
                    dbLog.addLogProp("... changed gene type in rgd", "EG_IN_RGD_DIFF_TYPE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
                }
            }
        }
        else {
            logger.debug("eg id and rgd id match, at least one sequence match, incoming eg type is not pseudo");
            dbLog.addLogProp("... incoming eg type is not pseudo", "UPDATE", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
        }
        return flag;
    }

    private boolean oneXdbMatch(List<XdbId> rgdXdbs, List<XdbId> newXdbs ) {

        for (XdbId rgdXdb: rgdXdbs) {
            for (XdbId newXdb: newXdbs) {
                if (newXdb.getAccId()!=null && rgdXdb.getAccId()!=null && newXdb.getAccId().trim().equals(rgdXdb.getAccId().trim())) {
                    return true;
                }
            }
        }        
        return false;        
    }
    
    
    public Flags orthologsCheck (BulkGene fileGene, XdbId xdbId) throws Exception {

        Flags flag;

        // check if the eg type is one of the types we want
        checkAllowedTypes(fileGene);

        EGDAO dao = EGDAO.getInstance();
        List<XdbId> nMHXdb = fileGene.getXdbIdsByXdbKey(xdbId.getXdbKey()); // the file MGD ID or HGNC ID
        List<Gene> genesByEgId = dao.getGenesByEGID(fileGene.getEgId()); //rgd genes matched by eg ID
        List newSeqXdbs = fileGene.getXdbIdsByXdbKey(XdbId.XDB_KEY_GENEBANKNU); // the GB sequence in the new record
        List<Gene> activeRgdGene = new ArrayList<>(); // will hold active rgd gene matched by entrezgene ID
        String info="";
        // get matched active gene by entrezgene id
        for (Gene rgdGeneTmp: genesByEgId) {
            // get the rgd gene status, and put the active gene into a list
            String rgdGeneStatus = dao.getRgdId(rgdGeneTmp.getRgdId()).getObjectStatus();
            if (rgdGeneStatus.equals("ACTIVE")) {
                activeRgdGene.add(rgdGeneTmp);                
                info += rgdGeneTmp.getRgdId()+",";
            }           
        }
        
        if (activeRgdGene.size()>1) {
            // Entrezgene ID matches more than one active entrezgene in RGD
            flag = new Flags(Flags.MULTIGENES, Flags.ERROR);
            flag.setRelatedInfo("Associated with multiple genes in RGD: "+ info);   
            dbLog.addLogProp(flag.getRelatedInfo(), "MULTIGENES", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }  
        else if (activeRgdGene.size()== 1){
            logger.debug("incoming EG is associated with one gene in RGD");
            flag = new Flags(Flags.EGINRGD);
            fileGene.rgdGene = activeRgdGene.get(0);
            dbLog.addLogProp("incoming EG is associated with one gene in RGD", "EG_IN_RGD", fileGene.getRecNo(), PipelineLogger.REC_FLAG);

            // check MGI id
            if (nMHXdb.size() >0) {
                // the new Entrezgene has the mgd id, check if two mgd id match
                //logger.debug("The eg id is in RGD, the incoming eg record has a MGD(HGNC) ID");
                String nMHAccId=nMHXdb.get(0).getAccId(); // the new MGD(HGNC) ID
                List<XdbId> rgdMHXdbList = dao.getXdbIdsByRgdId(xdbId.getXdbKey(), fileGene.rgdGene.getRgdId());
                boolean mhIdMatch=false;
                String rgdMgdId=null;
                // if only one MGD(HGNC) ID match, means the MGD(HGNC) ID match
                if (rgdMHXdbList.size() >0) {
                    for (XdbId accId: rgdMHXdbList) {
                        String rgdMgdIdTmp = accId.getAccId();

                        if (nMHAccId.trim().equals(rgdMgdIdTmp.trim())) {
                            mhIdMatch = true;
                            rgdMgdId = rgdMgdIdTmp;
                            break;
                        }
                        rgdMgdId = rgdMgdIdTmp;
                    }
                    if (!mhIdMatch) {
                        // the rgd mgd id and new mgd id doesn't match
                        logger.debug("The eg id is in RGD, incoming eg record has a MGD(HGNC) ID, two ID doesn't match");
                        flag = new Flags(Flags.EGINRGD_DIFFMHID, Flags.ERROR);
                        flag.setRgdId(fileGene.rgdGene.getRgdId());
                        flag.setRelatedInfo("RGD M(H) ID:"+rgdMgdId+" NEW M(H) ID:"+nMHAccId);
                        dbLog.addLogProp("eg id is in RGD, incoming eg record has a MGD(HGNC) ID, two ID doesn't match; "+flag.getRelatedInfo(), "EGINRGD_DIFFMHID", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                        return flag;
                    } else {
                        // log the matching MGD-ID or HGNC ID
                        dbLog.addLogProp(xdbId.getXdbKey()==XdbId.XDB_KEY_MGD?"MGD":"HGNC", rgdMgdId, fileGene.getRecNo(), PipelineLogger.REC_XDBID);
                        // continue to check sequence
                    }
                }
                else {
                    // the rgd gene doesn't have MGD ID
                    //check the sequence
                }
            }
            else {
                // the new entrezgene has no mgd id
                //check sequence
            }
        }
        else if( fileGene.getGene().getSpeciesTypeKey()>=4 ) {

            flag = new Flags(Flags.NEWGENE, Flags.INSERT);
            dbLog.addLogProp("new gene", "NEW_GENE", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }
        else {
            // the egid doesn't match active egid in rgd, start checking MGD (HGNC) ID             
            if (nMHXdb.size() >0) {
                // the new Entrezgene has the mgd id
                logger.debug("The eg id is not in RGD, incoming eg record has a MGD(HGNC) ID");
                String nMHAccId = nMHXdb.get(0).getAccId();
                List<Gene> genesByMHAccId = dao.getActiveGenesByXdbId(xdbId.getXdbKey(), nMHAccId); // get the gene by mgd or HGNC ID
                
                if (genesByMHAccId.size() ==1) {
                    // the new Entrezgene matches exactly one MGD ID in RGD
                    flag = new Flags(Flags.EGINRGD);
                    fileGene.rgdGene = genesByMHAccId.get(0); // will continue to check the sequence
                    dbLog.addLogProp("eg id is not in RGD, incoming eg record has a MGD(HGNC) ID", "EG_IN_RGD", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                    // log the matching MGD-ID or HGNC ID
                    dbLog.addLogProp(xdbId.getXdbKey()==XdbId.XDB_KEY_MGD?"MGD":"HGNC", nMHAccId, fileGene.getRecNo(), PipelineLogger.REC_XDBID);
                }
                else if (genesByMHAccId.size() >1) {
                    flag = new Flags(Flags.MULTIGENES, Flags.ERROR);
                    String relatedInfo="";
                    for (Gene rgdGeneTmp: genesByMHAccId) {
                        relatedInfo += "RGD:"+rgdGeneTmp.getRgdId()+",";
                    }
                    flag.setRelatedInfo(relatedInfo);

                    // log all the matching MGD-IDs or HGNC IDs
                    for( XdbId multiXdbId: nMHXdb ) {
                        dbLog.addLogProp(xdbId.getXdbKey()==XdbId.XDB_KEY_MGD?"MGD":"HGNC", multiXdbId.getAccId(), fileGene.getRecNo(), PipelineLogger.REC_XDBID);
                    }

                    dbLog.addLogProp("eg id is not in RGD, incoming eg record has multiple MGD(HGNC) ID: "+flag.getRelatedInfo(), "MULTIGENES", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                    return flag;
                }
                else {
                    // other species, eg id doesn't match, mgi id doesn't match, check if it has any sequence
                    if (newSeqXdbs.size() >0) {
                        // it has sequence, load as a new record
                        flag = new Flags(Flags.NEWGENE, Flags.INSERT);
                        dbLog.addLogProp("eg id doesn't match, mgi/hgnc id doesn't match, new sequence found", "NEW_GENE", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                        return flag;
                    } else {
                        // it has no sequence, go to no sequence log file
                        flag = new Flags(Flags.NOSEQ, Flags.ERROR);
                        dbLog.addLogProp("eg id doesn't match, mgi/hgnc id doesn't match, no sequence found", "NO_SEQ", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                        return flag;
                    }
                }
            }
            else {
                // new gene has no mgd or HGNC id
                flag = new Flags(Flags.NOMHID, Flags.SKIP);
                dbLog.addLogProp("new record doesn't have mgi/hgnc ID", "NO_MHID", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                return flag;
            }
        }        
        
        /*
         to here, entrezgene ID match active gene in RGD or MGIID match active gene in RGD        
         the flag value is EGINRGD     
        */
        
        //check if the MGD id match        
        flag.setRgdId(fileGene.rgdGene.getRgdId());
        flag.setRelatedInfo("");  
          
        /*
         * Check the sequence
         */
        List<XdbId> rgdSeqXdbs = dao.getXdbIdsByRgdId(XdbId.XDB_KEY_GENEBANKNU, flag.getRgdId());
        rgdSeqXdbs = xdbManager.removeGeneBankNucleotides(rgdSeqXdbs);
        if (newSeqXdbs.size() ==0) {
            // exception: do not enforce sequence requirement for mitochondrial genes
            if( !Utils.stringsAreEqual(fileGene.getChromosome(),"MT") ) {
                flag = new Flags(Flags.NOSEQ, Flags.ERROR, flag.getRgdId());
                flag.setRelatedInfo("RGD:"+ fileGene.rgdGene.getRgdId());
                dbLog.addLogProp("no seq"+flag.getRelatedInfo(), "NO_SEQ", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                return flag;
            }
        }

        if (rgdSeqXdbs.size() <=0 ) {
            // the rgd record has no sequence
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp("rgd record has no sequence", "UPDATE", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }

        if (newSeqXdbs.size() <= 0) {
            flag = new Flags(Flags.NOSEQ, Flags.SKIP, flag.getRgdId());
            dbLog.addLogProp("no sequence", "NO_SEQ", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }

        //the new entrezgene record has at least one sequence
        if (oneXdbMatch(rgdSeqXdbs, newSeqXdbs)) {
            // the eg id and rgd id match, at least one sequence match proceed to update
            flag.setLoadStatus(Flags.UPDATE);
            dbLog.addLogProp("the eg id and rgd id match, at least one sequence match proceed to update", "UPDATE", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
            return flag;
        }
        else {
            // is there a single XM or XR sequence in RGD?
            boolean isXmXr=false;
            for (XdbId rgdXdb: rgdSeqXdbs) {
                if (rgdXdb.getAccId()!=null) {
                    if (rgdXdb.getAccId().startsWith("XM_") || rgdXdb.getAccId().startsWith("XR_")){
                        isXmXr=true;
                    }
                }
            }
            // there is XM  or XR sequence in RGD
            if (isXmXr){
            // the eg id and rgd id match, sequence doesn't match, has XM sequence
                logger.debug("the eg id and rgd id match, sequence doesn't match, has XM sequence");
                flag.setLoadStatus(Flags.UPDATE);
                dbLog.addLogProp("the eg id and rgd id match, sequence doesn't match, has XM sequence", "UPDATE", fileGene.getRecNo(), PipelineLogger.REC_FLAG);
                return flag;
            }
            else {
                // the eg id and rgd id match, sequence doesn't match, no sequence is XM or XR
                handleMissingSequence(flag, fileGene, rgdSeqXdbs, newSeqXdbs);
                return flag;
            }
        }
    }
     
    private String rgdAccXdbToString(List<XdbId> rgdAccXdb) throws Exception {
        return Utils.concatenate(",", rgdAccXdb, "getAccId");
    }

    Flags checkEgId(BulkGene bulkGene ) {

        // log the RGDID as read from EG record
        dbLog.addLogProp("found RGDID", Integer.toString(bulkGene.getGene().getRgdId()), bulkGene.getRecNo(), PipelineLogger.REC_RGDID);

        // check if the new EntrezGene has an EntrezGene ID
        if (!bulkGene.hasEgId() ) {
            // new bulk gene doesn't have an entrezgene ID
            dbLog.addLogProp("no EntrezGene ID - record skipped", "NO_EG_ID", bulkGene.getRecNo(), PipelineLogger.REC_FLAG);
            return new Flags(Flags.NOEGID, Flags.ERROR);
        }

        logger.debug("incoming EG:" + bulkGene.toString());
        dbLog.addLogProp("EGID", bulkGene.getEgId(), bulkGene.getRecNo(), PipelineLogger.REC_XDBID);
        return null;
    }

    void checkAllowedTypes(BulkGene bg ) throws Exception {

        // check if incoming gene type is among the allowed gene types
        if (!excludeEgType.contains(bg.getEgType())) {
            return;
        }

        // not allowed gene type: flag the record and change gene type to 'gene'
        bg.getSession().incrementCounter("WRONG_GENE_TYPE", 1);
        getDbFlagManager().setFlag("WRONG_GENE_TYPE", bg.getRecNo());
        bg.setEgType("gene");
    }

    // STAGE 2: additional check after the standard quality check had been run
    public void afterCheck (BulkGene bg) throws Exception {

        Flags flag = bg.getCustomFlags();

        // both the eg full name and rgd full name must be present
        // to conduct full name comparison
        String egGeneName = bg.getGene().getName();
        String rgdGeneName = bg.rgdGene!=null ? bg.rgdGene.getName() : "";
        if( !Utils.isStringEmpty(egGeneName) && !Utils.isStringEmpty(rgdGeneName)) {
            if( !egGeneName.equals(rgdGeneName.trim())) {
                String msg = "eg gene name differs from rgd full name: EG-NAME="+egGeneName+", RGD-NAME:"+rgdGeneName;
                dbLog.addLogProp(msg, "GENE_NAME_MISMATCH", bg.getRecNo(), PipelineLogger.REC_FLAG);

                // modify flags
                flag.setFlagValue(flag.getFlagValue()+","+ Flags.GENENAME_MISMATCH);
                flag.setRelatedInfo(flag.getRelatedInfo()+","+msg);
            }
        }

        // compute status of transcripts BEFORE exporting bulkgene to XML
        computeCodingStatusForTranscripts(bg);

        if(flag.getLoadStatus().equals(Flags.UPDATE)) {
            // get a copy of rgd map data now, so data load thread will have the data ready
            // note: include only those map entries that are on official list of maps processed by pipeline
            bg.genePositions.loadRgdMapData(flag.getRgdId(), validMapKeys);

            bg.genePositions.qcMapData(bg, null);
        }
        else if( flag.getLoadStatus().equals(Flags.INSERT) ) {
            bg.genePositions.qcMapData(bg, null);
        }

        // check for changes in gene symbol, gene name and gene description
        checkForNomenclatureChanges(bg.getGene(), bg.rgdGene, bg);

        // recreate xml representation
        bg.setXmlString(bg.toXmlString());

        // check if gene refseq status changed
        checkForRefSeqChanges(bg);

        // check if gene ncbi annot status changed
        checkForNcbiAnnotStatusChanges(bg);

        // gather data for data loading
        if(flag.getLoadStatus().equals(Flags.UPDATE)) {

            int rgdId = flag.getRgdId();

            bg.aliases.qcIncomingAliases(bg.rgdGene, rgdId, bg.getSession());

            qcXdbIds(bg, rgdId);

            //##### TRANSCRIPTS
            // load rgd transcripts
            bg.rgdTranscripts = transcriptDAO.getTranscriptsForGene(rgdId);
            removeTranscriptsWithInvalidMapKeys(bg.rgdTranscripts);
            // load all genomic features with positions
            bg.rgdFeatures = transcriptDAO.getFeaturesForGene(rgdId);
            removeFeaturesWithInvalidMapKeys(bg.rgdFeatures);
        }
        else
        if(flag.getLoadStatus().equals(Flags.INSERT)) {
            // all aliases need to be added
            bg.aliases.qcIncomingAliases(null, 0, bg.getSession());
        }

        // none of the aliases can have a value identical to gene name or symbol
        bg.aliases.qcAliases(bg);

        // initialize transcript data, if not done yet
        if( bg.rgdTranscripts==null )
            bg.rgdTranscripts = new ArrayList<Transcript>();
        if( bg.rgdFeatures==null )
            bg.rgdFeatures = new ArrayList<TranscriptFeature>(32);
    }

    void qcXdbIds(BulkGene bg, int rgdId) throws Exception {

        // get a copy of rgd xdb ids for given rgdid and ENTREZGENE pipeline
        XdbId xdbId = new XdbId();
        xdbId.setRgdId(rgdId);
        xdbId.setSrcPipeline(XdbManager.EG_PIPELINE);
        List<XdbId> rgdXdbIds = bg.dao.getXdbIds(xdbId);

        // set rgdid for incoming xdb ids
        List<XdbId> xdbs = bg.getXdbIds();
        if (xdbs!=null && xdbs.size() >0 ) {
            // ensure all incoming xdb ids have valid rgd-id assigned
            for( XdbId xdb: xdbs ) {
                xdb.setRgdId(rgdId);
            }
        }
        // construct a set of duplicates: xdbids that are present both in rgd and in incoming list
        List<XdbId> matchingXdbs = new ArrayList<XdbId>(rgdXdbIds);
        matchingXdbs.retainAll(xdbs);
        // remove from incoming xdbids those that are already in rgd
        List<XdbId> xdbsIncoming = new ArrayList<XdbId>(xdbs);
        xdbsIncoming.removeAll(rgdXdbIds);
        // remove from rgdids the duplicates (present on incoming xdbid list)
        rgdXdbIds.removeAll(matchingXdbs);
        bg.getSession().incrementCounter("XDBIDS_MATCHED", matchingXdbs.size());
        // update bg
        bg.toBeInsertedXdbIds = xdbsIncoming;
        bg.toBeRemovedXdbIds = rgdXdbIds;
        bg.matchingXdbIds = matchingXdbs;

        qcKeggPathwayIds(bg);
    }

    void qcKeggPathwayIds(BulkGene bg) throws Exception {

        // if there are any KEGG pathway ids among to-be-inserted or to-be-removed array,
        // they should be removed
        Iterator<XdbId> it = bg.toBeInsertedXdbIds.iterator();
        while( it.hasNext() ) {
            if( it.next().getXdbKey()==XdbId.XDB_KEY_KEGGPATHWAY )
                it.remove();
        }

        it = bg.toBeRemovedXdbIds.iterator();
        while( it.hasNext() ) {
            if( it.next().getXdbKey()==XdbId.XDB_KEY_KEGGPATHWAY )
                it.remove();
        }

        KeggPathwayIdsSyncer syncer = new KeggPathwayIdsSyncer();
        syncer.dao = bg.dao;

        // get incoming KEGG pathway ids
        List<XdbId> incomingKeggPathwayIds = bg.getXdbIdsByXdbKey(XdbId.XDB_KEY_KEGGPATHWAY);
        for( XdbId xdbId: incomingKeggPathwayIds ) {
            syncer.addIncomingObject(xdbId);
        }

        syncer.qc(bg.getCustomFlags().getRgdId());

        bg.toBeInsertedXdbIds.addAll(syncer.getForInsertList());
        bg.toBeUpdatedXdbIds = syncer.getForUpdateList();

        // never remove ids not created by 'ENTREZGENE' pipeline
        for( XdbId xdbId: (List<XdbId>) syncer.getForDeleteList() ) {
            if( Utils.stringsAreEqualIgnoreCase(xdbId.getSrcPipeline(), XdbManager.EG_PIPELINE) ) {
                bg.toBeRemovedXdbIds.add(xdbId);
            }
        }
    }

    void computeCodingStatusForTranscripts(BulkGene bg) {

        for( TranscriptInfo ti: bg.transcripts ) {
            ti.setNonCoding(ti.getProteinAccId()==null || ti.getProteinAccId().isEmpty() );
        }
    }

    void checkForRefSeqChanges(BulkGene bg) {

        // check if gene refseq status changed
        if( bg.rgdGene!=null && !Utils.stringsAreEqualIgnoreCase(bg.getGene().getRefSeqStatus(), bg.rgdGene.getRefSeqStatus())) {
            String msg = "gene refseq status changed: OLD="+bg.rgdGene.getRefSeqStatus()+", NEW:"+bg.getGene().getRefSeqStatus();
            dbLog.addLogProp(msg, "GENE_REFSEQ_CHANGED", bg.getRecNo(), PipelineLogger.REC_FLAG);

            // modify flags
            bg.setFlag("GENE_REFSEQ_STATUS");
        }
    }

    void checkForNcbiAnnotStatusChanges(BulkGene bg) {

         Flags flag = bg.getCustomFlags();

         // check if gene ncbi annot status changed
         if( bg.rgdGene!=null && !Utils.stringsAreEqualIgnoreCase(bg.getGene().getNcbiAnnotStatus(), bg.rgdGene.getNcbiAnnotStatus())) {
             String msg = "gene ncbi annotation status changed: OLD={"+bg.rgdGene.getNcbiAnnotStatus()+"}, NEW:{"+bg.getGene().getNcbiAnnotStatus()+"}";
             dbLog.addLogProp(msg, "NCBI_ANNOT_STATUS", bg.getRecNo(), PipelineLogger.REC_FLAG);

             // modify flags
             flag.setFlagValue(flag.getFlagValue()+","+ Flags.NCBI_ANNOT_STATUS);
             flag.setRelatedInfo(flag.getRelatedInfo()+","+msg);
         }
     }

    // this should be called only for orthologs: mouse and human genes
    void checkForNomenclatureChanges(Gene egGene, Gene rgdGene, BulkGene bg) {

        if( rgdGene==null )
            return;

        if( !Utils.stringsAreEqual(egGene.getSymbol(), rgdGene.getSymbol())) {
            logNomen.info("gene symbol change: RGD_ID=[" + rgdGene.getRgdId() + "], OLD=[" + rgdGene.getSymbol() + "], NEW=[" + egGene.getSymbol() + "]");
            logSymbol.info("|"+SpeciesType.getCommonName(rgdGene.getSpeciesTypeKey())
                    +"|RGDID:"+rgdGene.getRgdId()
                    +"|RGD SYMBOL|"+rgdGene.getSymbol()
                    +"|NCBI SYMBOL|"+egGene.getSymbol());
            bg.setFlag("ORTHO_NOMEN_SYMBOL");

            Alias alias = new Alias();
            alias.setValue(rgdGene.getSymbol());
            alias.setTypeName("old_gene_symbol");
            alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
            bg.aliases.addAlias(alias);
        }

        if( !Utils.stringsAreEqual(egGene.getName(), rgdGene.getName())) {
            logNomen.info("gene name change: RGD_ID=["+rgdGene.getRgdId()+"], OLD=["+rgdGene.getName()+"], NEW=["+egGene.getName()+"]");
            bg.setFlag("ORTHO_NOMEN_NAME");

            Alias alias = new Alias();
            alias.setValue(rgdGene.getName());
            alias.setTypeName("old_gene_name");
            alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
            bg.aliases.addAlias(alias);
        }

        if( !Utils.stringsAreEqual(egGene.getDescription(), rgdGene.getDescription())) {
            logNomen.info("gene description change: RGD_ID=["+rgdGene.getRgdId()+"], OLD=["+rgdGene.getDescription()+"], NEW=["+egGene.getDescription()+"]");
            bg.setFlag("ORTHO_NOMEN_DESC");
        }
    }

    // remove transcripts having positions on assembly maps not processed by the pipeline
    void removeTranscriptsWithInvalidMapKeys(List<Transcript> trs) {
        for (Transcript tr : trs) {
            // remove all genomic positions that are not a list of maps processed by the pipeline
            Iterator<MapData> mdIt = tr.getGenomicPositions().iterator();
            while (mdIt.hasNext()) {
                MapData md = mdIt.next();
                if (!validMapKeys.contains(md.getMapKey())) {
                    mdIt.remove();
                }
            }
        }
    }

    // remove transcript features having positions on assembly maps not processed by the pipeline
    void removeFeaturesWithInvalidMapKeys(List<TranscriptFeature> trfs) {
        Iterator<TranscriptFeature> it = trfs.iterator();
        while( it.hasNext() ) {
            TranscriptFeature tf = it.next();

            // remove the feature if it s located on assembly map not processed by the pipeline
            if( !validMapKeys.contains(tf.getMapKey()) )
                it.remove();
        }
    }

    public List getEgType() {
        return egType;
    }
    public void setEgType(List egType) {
        this.egType = egType;
    }    
    
    public List getExcludeEgType() {
        return excludeEgType;
    }
    public void setExcludeEgType(List excludeEgType) {
        this.excludeEgType = excludeEgType;
    }

    public Map<String, String> getGenomicAssemblies() {
        return genomicAssemblies;
    }

    public void setGenomicAssemblies(Map<String, String> genomicAssemblies) {
        this.genomicAssemblies = genomicAssemblies;

        // after genomic assemblies has been set, the valid map keys must be calculated
        validMapKeys = new HashSet<>();
        if( genomicAssemblies!=null ) {
            for (String val : genomicAssemblies.values()) {
                validMapKeys.add(Integer.parseInt(val));
            }
        }
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public PipelineLogFlagManager getDbFlagManager() {
        return dbFlagManager;
    }

    public void setDbFlagManager(PipelineLogFlagManager dbFlagManager) {
        this.dbFlagManager = dbFlagManager;
    }
}
