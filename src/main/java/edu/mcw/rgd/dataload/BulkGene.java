package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.*;
import edu.mcw.rgd.xml.BulkGeneXmlCreator;

import java.util.*;
import java.util.Map;

/**
 * contains all information extracted by analysis of a single EntrezGene record
 */
public class BulkGene {
    int recNo;
	Gene gene = new Gene();
    String egType;
    String chromosome, fishband;
    Double cM;
    private String egId; // EntrezGene ID
	Flags flags = new Flags();
    String xmlString;
    EGDAO dao;
    Set<String> _flags = new HashSet<>();

    // gene track status - as of September 2010, a gene can be "live", "replaced", "discontinued" or "unknown"
    // LIVE - The record is current and primary, i.e., not secondary or discontinued.
    // SECONDARY - The record is no longer current because it has been made secondary to another gene record.
    //     (The term secondary is applied to any record that has been merged into another.
    //      This occurs most often when multiple genes are defined based on incomplete data, and these are later
    //      discovered to be parts of the same gene. One gene record then becomes secondary to the other.)
    // DISCONTINUED - The record is no longer current, and it has not been made secondary to any other gene record.
    // UNKNOWN - No data available about gene track status.
    String geneTrackStatus;
    String geneTrackCurrentId; // primary gene id for a secondary record

    public String chromosome2, bioSourceGenome;
    public String nomenSymbol, nomenName; // gene official nomenclature name and symbol
    public String geneOfficialSymbol, geneOfficialName;
    public String geneInterimSymbol, geneInterimName;
    // NOTE: nomenSymbol/Name, geneOfficialSymbol/Name, geneInterimSymbol/Name are auxiliary fields
    //  to set Name and Symbol properties of gene object when name and symbol could not be determined
    // from incoming data; after XML parsing stage, gene.Name and .Symbol properties represent the incoming data

    // only set for regulatory elements
    public String biologicalRegionType;

    public Gene rgdGene; // set during quality check when EG ID could be connected to an active rgd-Gene
    public GenomicElement rgdElement; // set during quality check when EG ID could be connected to an active genomic element in RGD

    // handles the aliases
    public AliasLoader aliases = new AliasLoader();

    // auxiliary data evaluated during second phase of quality check
    // (if evaluated here, the data loading will be faster)

    private List<XdbId> xdbIds = new ArrayList<>(); // incoming xdbs

    public GenePositions genePositions;

    // transcript data in rgd;
    // rgd transcripts
    public List<Transcript> rgdTranscripts;
    // all genomic features with positions in rgd
    public List<TranscriptFeature> rgdFeatures;
    //
    public TranscriptXdbIds transcriptXdbIds;

    // temporary data valid only during xml parsing:
    // map of transcript accession id to refseq_status
    public Map<String, String> mapRefSeqStatus = new HashMap<>();

	public BulkGene(){
		gene.setSymbol("empty");
        dao = EGDAO.getInstance();
        genePositions = new GenePositions(dao);
        transcriptXdbIds = new TranscriptXdbIds();
	}

    // map key for primary reference assembly; should be set at the very beginning of processing
    static int _primaryMapKey = SpeciesType.ALL;
    static public int getPrimaryMapKey() {
        return _primaryMapKey;
    }
    static public void setPrimaryMapKey(int mapKey) {
        _primaryMapKey = mapKey;
    }

    public EGDAO getDao() {
        return dao;
    }

	public Flags getCustomFlags() {
		return flags;
	}
    
	public void setCustomFlags(Flags flags) {
		this.flags = flags;
	}
	
	public Gene getGene() {
		return gene;
	}
	
    public List<XdbId> getXdbIds() {
        return this.xdbIds;
    }

    // add new record only if not a duplicate
	public void addXdbId(XdbId newRec){

        for( XdbId oldRec: xdbIds ) {
            if( oldRec.equals(newRec) )
                return; // newRec already exists -- we won't add it twice
        }
		this.xdbIds.add(newRec);
	}

	public List<XdbId> getXdbIdsByXdbKey(int xdbKey) {

		ArrayList<XdbId> xdbs = new ArrayList<> ();
        for (XdbId xdbId: xdbIds) {
            if( xdbId.getXdbKey()==xdbKey ) {
                if (xdbId.getAccId()!=null)
                    xdbs.add(xdbId);
            }
        }
        return xdbs;
    }

    public void removeObsoleteHgncIds(EGDAO dao) {

        xdbIds.removeIf( id -> {
            try {
                return id.getXdbKey()==XdbId.XDB_KEY_HGNC && dao.isObsoleteHgncId(id.getAccId());
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    public boolean hasRefSeqNucleotide() {

        for( XdbId xdbId: getXdbIdsByXdbKey(XdbId.XDB_KEY_GENEBANKNU) ) {
            if( xdbId.getAccId().startsWith("NM_") ||
                    xdbId.getAccId().startsWith("XM_") ||
                    xdbId.getAccId().startsWith("NG_") ||
                    xdbId.getAccId().startsWith("NR_") ||
                    xdbId.getAccId().startsWith("XR_") ) {
                return true;
            }
        }
        return false;
    }

    public boolean hasRefSeqCuratedStatus() {

        return Utils.stringsAreEqualIgnoreCase(this.gene.getRefSeqStatus(), "REVIEWED") ||
                Utils.stringsAreEqualIgnoreCase(this.gene.getRefSeqStatus(), "VALIDATED");
    }

	public String toString(){
        return this.getEgId() + '\t' + this.getEgType() +
                '\t' + this.getGene().getRgdId() +
                '\t' + this.getGene().getSymbol();
	}
    
    public String toXmlString() {
        BulkGeneXmlCreator bgXmlCreator = BulkGeneXmlCreator.getInstance();
        return bgXmlCreator.createXmlString(this);
    }
    
    public String getEgType() {
        return egType;
    }

    public void setEgType(String egType) {
        this.egType = egType;
    }

    // return EG-ID as string (always non-null)
    public String getEgId() {
        return egId==null ? "" : egId;
    }
    public void setEgId(String egId) {
        this.egId = egId;
    }

    // true if there is a valid Entrez Gene ID
    public boolean hasEgId() {
        String id = getEgId();
        if( Utils.isStringEmpty(id) )
            return false;
        return Integer.parseInt(getEgId()) > 0;
    }

    public String getXmlString() {
        return xmlString;
    }
    public void setXmlString(String xmlString) {
        this.xmlString = xmlString;
    }
    
    public String getChromosome() {
        return chromosome;
    }
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }
    
    public String getFishband() {
        return fishband;
    }
    public void setFishband(String fishband) {
        this.fishband = fishband;
    }

    public Double getcM() {
        return cM;
    }
    
    public void setcM(Double cM) {
        this.cM = cM;
    }

    /**
     * represents a transcript and a set of genomic features, like exons, utrs, introns, and so on
     */
    public LinkedList<TranscriptInfo> transcripts = new LinkedList<TranscriptInfo>();

    public TranscriptInfo getTranscriptByAccId(String accId) {

        for( TranscriptInfo tr: transcripts ) {
            if( Utils.stringsAreEqual(tr.getAccId(), accId) ) {
                return tr;
            }
        }
        return null;
    }

    public TranscriptInfo createTranscript(String accId) {

        if( Utils.isStringEmpty(accId) ) {
            // attempt to add a null transcript
            System.out.println("attempt to add a null transcript");
            return null;
        }

        TranscriptInfo tr = new TranscriptInfo();
        tr.setAccId(accId);
        transcripts.add(tr);
        return tr;
    }

    public void addMapDataForGene(MapData mdNew) {

        // add new position
        genePositions.addMapData(mdNew);
    }

    public boolean hasGenomicPositionForAssembly(int mapKey) {
        return genePositions.hasGenomicPositionForAssembly(mapKey);
    }

    public void setNcbiAnnotStatus(String ncbiAnnotStatus) {
        // multiple annotation statuses could be present for one gene; combine them all
        String newStatus = this.gene.getNcbiAnnotStatus();
        if( newStatus==null )
            newStatus = ncbiAnnotStatus;
        else
            newStatus += "; " + ncbiAnnotStatus;
        this.gene.setNcbiAnnotStatus(newStatus);
    }

    public String getGeneTrackStatus() {
        return geneTrackStatus;
    }

    public void setGeneTrackStatus(String geneTrackStatus) {
        this.geneTrackStatus = geneTrackStatus;
    }

    public String getGeneTrackCurrentId() {
        return geneTrackCurrentId;
    }

    public void setGeneTrackCurrentId(String geneTrackCurrentId) {
        this.geneTrackCurrentId = geneTrackCurrentId;
    }

    public int getRecNo() {
        return recNo;
    }

    public void setRecNo(int recNo) {
        this.recNo = recNo;
    }

    public void setFlag(String flag) {
        _flags.add(flag);
    }

    public boolean isFlagSet(String flag) {
        return _flags.contains(flag);
    }
}


