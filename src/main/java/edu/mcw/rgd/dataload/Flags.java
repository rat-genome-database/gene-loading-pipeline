package edu.mcw.rgd.dataload;

// flags corresponding to different stages of analysis for the current pipeline
public class Flags {
    public static final String NEWGENE= "NEW_GENE";
    public static final String NONACTIVE= "NON_ACTIVE_IN_RGD";
    public static final String DIFFEGID="EGID_NOT_IN_RGD";   
    public static final String EGINRGD="EGID_IN_RGD";
    public static final String EGINRGD_DIFFRGDID="EG_IN_RGD_DIFF_RGDID";
    public static final String EGINRGD_DIFFMHID="EG_IN_RGD_DIFF_MGDID(HGNCID)";
    public static final String EGINRGD_VAR="EGID_IN_RGD_VARIATION";

    public static final String NOMHID="NO_M-H_ID";
    public static final String NOEGID="NO_ENTREZGENE_ID";    
    public static final String MULTIGENES="MULTIPLE_GENES";
    public static final String GENENAME_MISMATCH="GENE_NAME_MISMATCH";
    public static final String NCBI_ANNOT_STATUS="NCBI_ANNOT_STATUS";

    public static final String INSERT="INSERT";
    public static final String UPDATE="UPDATE";
    public static final String ERROR="CONFLICT"; 
    public static final String SKIP="SKIP";

	int rgdId;
    String flagValue;
    String loadStatus;
	String relatedInfo;    
	
    public Flags() {
    }

	// create a flag with a given flagValue
    public Flags(String flagValue) {
        this.flagValue = flagValue;
    }

    // create a flag with a given flagValue and a loadStatus
    public Flags(String flagValue, String loadStatus) {
        this.flagValue = flagValue;
        this.loadStatus = loadStatus;
    }

    public String getFlagValue() {
		return flagValue;
	}
	public void setFlagValue(String flagValue) {
		this.flagValue = flagValue;
	}
	public String getRelatedInfo() {
		return relatedInfo;
	}
	public void setRelatedInfo(String relatedInfo) {
		this.relatedInfo = relatedInfo;
	}
    
    public String getLoadStatus() {
        return loadStatus;
    }
    public void setLoadStatus(String loadStatus) {
        this.loadStatus = loadStatus;
    }
    
    public int getRgdId() {
        return rgdId;
    }
    public void setRgdId(int rgdId) {
        this.rgdId = rgdId;
    }
    public String toString() {
        return rgdId+"\t" + flagValue.concat("\t"+ loadStatus).concat("\t"+ relatedInfo);
    }
}

