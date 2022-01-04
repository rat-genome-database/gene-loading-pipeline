package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.PipelineLogDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogFlagManager;
import edu.mcw.rgd.process.PipelineLogger;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Date;

/**
 * Read the flag in the Bulk Gene, and decide if it is insert, or update or write to log file
 */
public class DecisionMaker {
    
    BulkGeneLoaderImpl bulkGeneLoader;
    int newGenes=0;
    int updated=0;
    int error=0;
    int noMHID=0;
    int skipped=0;
    int nameMismatch=0;
    int geneRefSeq=0; // count of genes with changed refseq status
    PipelineLogger dbLog = null;
    PipelineLogFlagManager dbFlagManager;

    protected final Logger logger = LogManager.getLogger("process");

    // look at the flags and decide if the record need update or insert
    public void decide(BulkGene bg, CounterPool counters) throws Exception {

        bulkGeneLoader.setCounters(counters);

        dbLog = PipelineLogger.getInstance();

        if( bg.isFlagSet("BAD_SPECIES") ) {
            dbFlagManager.setFlag("BAD_SPECIES", bg.getRecNo());
            decide(bg, bg.getCustomFlags());
        } else {
            decide(bg, bg.getCustomFlags());
            handleGeneTrackStatus(bg);
        }

        dbFlagManager.writeFlags(bg.getRecNo());
    }

    public void decide(BulkGene bg, Flags flag) throws Exception {

        // write XML representation of the record to log file
        String xml = bg.getXmlString();
        if( xml!=null ) {
            String info = "";
            String value = "";

            PipelineLog plog = dbLog.getPipelineLog();
            PipelineLogDAO plogDAO = dbLog.getPipelineLogDAO();
            PipelineLog.LogProp xmlProp = plog.createLogProp();
            xmlProp.setKey(PipelineLogger.REC_XML);
            xmlProp.setDate(new Date());
            xmlProp.setInfo(info);
            xmlProp.setValue(value);
            xmlProp.setXml(xml);
            xmlProp.setRecNo(bg.getRecNo());
            plogDAO.writeLogProp(plog.getPipelineLogKey(), xmlProp);
        }

        // write to gene variation log and diff type log
        if (flag.getFlagValue()==null || flag.getLoadStatus()==null) {
            logger.error("Not flagged in Quality check");
            return;
        }

        if (flag.getFlagValue().contains(Flags.GENENAME_MISMATCH)) {
            nameMismatch++;
        }

        // update gene refseq status
        if (bg.isFlagSet("GENE_REFSEQ_STATUS")) {
            geneRefSeq++;
        }

        // insert, update or upgrade depends on the flag value
        if (flag.getFlagValue().equals(Flags.NOMHID)) {
            noMHID ++;
        } else if (flag.getLoadStatus().equals(Flags.SKIP)) {
            throw new Exception("Unexpected SKIP status");
        }


        // insert, update or upgrade depends on the flag value
        switch (flag.getLoadStatus()) {
            case Flags.ERROR:
                error++;
                break;
            case Flags.INSERT:
                // do not insert new genes that are flagged as DISCONTINUED or SECONDARY (very rare, but could happen)
                if (bg.getGeneTrackStatus().equals("DISCONTINUED") || bg.getGeneTrackStatus().equals("SECONDARY") ) {
                    skipped++;
                    dbFlagManager.setFlag("NEW_GENE_SKIPPED", bg.getRecNo());
                }
                // also do not insert genes with wrong type: 'other' or 'unknown'
                //   unless they have at least one RefSeq nucleotide: NM_, XM_, NR_, XR_, NG_
                else if( Utils.stringsAreEqualIgnoreCase(bg.getEgType(), "gene")
                        && !bg.hasRefSeqNucleotide()
                        && !bg.hasRefSeqCuratedStatus() ){
                    skipped++;
                    dbFlagManager.setFlag("NEW_GENE_SKIPPED", bg.getRecNo());
                } else {
                    bulkGeneLoader.load(bg);
                    newGenes++;
                    dbFlagManager.setFlag("NEW_GENE_INSERTED", bg.getRecNo());
                }
                break;
            case Flags.UPDATE:
                bulkGeneLoader.update(bg);
                updated++;
                break;
        }
    }

    void handleGeneTrackStatus(BulkGene bg) throws Exception {
        dbFlagManager.setFlag("GENE_TRACK_STATUS_"+bg.getGeneTrackStatus(), bg.getRecNo());

        // special logging for unusual gene track status
        if( !bg.getGeneTrackStatus().equals("LIVE") ) {
            int speciesTypeKey = bg.getGene()!=null ? bg.getGene().getSpeciesTypeKey() : DataLoadingManager.getInstance().getSpeciesTypeKey();
            String msg = "GeneTrackStatus="+bg.getGeneTrackStatus()
                    +"|OldGeneId="+bg.getEgId()
                    +"|GeneRGDId="+bg.getCustomFlags().getRgdId()
                    +(bg.getGene()!=null ? "|GeneSymbol="+bg.getGene().getSymbol() : "|")
                    +"|CurrentGeneId="+bg.getGeneTrackCurrentId()
                    +"|Species="+SpeciesType.getCommonName(speciesTypeKey);

            // UNKNOWN gene track status in the past years was associated entirely with mitochondrial genes
            // so we don't want to report this as problematic
            Logger log = LogManager.getLogger("geneTrackStatus");
            if( bg.getGeneTrackStatus().equals("UNKNOWN") ) {
                log.debug(msg);
            } else {
                log.info(msg);
            }

            // log unusual gene status with gene track status tracker
            DataLoadingManager.getInstance().getGeneStatusIssueTracker().addIssue(
                    bg.getGeneTrackStatus(),
                    bg.getCustomFlags().getRgdId(),
                    speciesTypeKey,
                    bg.getGene()!=null ? bg.getGene().getSymbol() : "",
                    bg.getEgId(),
                    bg.getGeneTrackCurrentId()
            );
        }
    }

    public int getError() {
        return error;
    }
    public int getNewGenes() {
        return newGenes;
    }
    public int getUpdated() {
        return updated;
    }
    public int getNoMHID() {
        return noMHID;
    }
    
    public BulkGeneLoaderImpl getBulkGeneLoader() {
        return bulkGeneLoader;
    }

    public void setBulkGeneLoader(BulkGeneLoaderImpl bulkGeneLoader) {
        this.bulkGeneLoader = bulkGeneLoader;
    }

    public int getNameMismatch() {
        return nameMismatch;
    }

    public void setNameMismatch(int nameMismatch) {
        this.nameMismatch = nameMismatch;
    }

    public PipelineLogFlagManager getDbFlagManager() {
        return dbFlagManager;
    }

    public void setDbFlagManager(PipelineLogFlagManager dbFlagManager) {
        this.dbFlagManager = dbFlagManager;
        if( this.bulkGeneLoader!=null )
            this.bulkGeneLoader.setDbFlagManager(dbFlagManager);
    }
}
