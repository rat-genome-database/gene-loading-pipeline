package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Set;
import java.util.concurrent.ConcurrentSkipListSet;

/**
 * Created by mtutaj on 2/2/2016.
 * <p>
 * gene track status issues, as reported by QC module, are being tracked here
 */
public class GeneStatusIssueTracker {

    private Set<GeneStatusIssue> geneStatusIssues = new ConcurrentSkipListSet<>();

    public void addIssue(String geneTrackStatus, int geneRgdId, int speciesTypeKey, String geneSymbol, String oldGeneId, String newGeneId) {
        GeneStatusIssue issue = new GeneStatusIssue();
        issue.geneTrackStatus = geneTrackStatus;
        issue.geneRgdId = geneRgdId;
        issue.speciesTypeKey = speciesTypeKey;
        issue.geneSymbol = geneSymbol;
        issue.oldGeneId = oldGeneId;
        issue.newGeneId = newGeneId;
        addIssue(issue);
    }

    public void addIssue(GeneStatusIssue issue) {
        geneStatusIssues.add(issue);
    }

    public void writeIssuesToFile() throws IOException {

        String fileName = "logs/" + SpeciesType.getCommonName(DataLoadingManager.getInstance().getSpeciesTypeKey());
        fileName += "GeneStatusIssues_";
        fileName += new SimpleDateFormat("yyyyMMdd").format(new Date());
        fileName += ".txt";

        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));

        // write header
        out.write("gene track status\tgene rgd id\tspecies type key\tgene symbol\told gene id\tnew gene id\n");

        // write all issues
        for( GeneStatusIssue issue: geneStatusIssues ) {
            out.write(issue.geneTrackStatus+'\t'+
                issue.geneRgdId+'\t'+
                issue.speciesTypeKey+'\t'+
                issue.geneSymbol+'\t'+
                issue.oldGeneId+'\t'+
                issue.newGeneId+'\n');
        }
        out.close();
    }

    public class GeneStatusIssue implements Comparable<GeneStatusIssue> {
        public String geneTrackStatus; // DISCONTINUED, SECONDARY
        public int geneRgdId;
        public int speciesTypeKey;
        public String geneSymbol;
        public String oldGeneId;
        public String newGeneId;

        @Override
        public int compareTo(GeneStatusIssue o) {
            int r = Utils.stringsCompareTo(this.geneTrackStatus, o.geneTrackStatus);
            if( r!=0 )
                return r;
            r = Utils.stringsCompareTo(this.geneSymbol, o.geneSymbol);
            if( r!=0 )
                return r;
            r = Utils.stringsCompareTo(this.oldGeneId, o.oldGeneId);
            if( r!=0 )
                return r;
            r = Utils.stringsCompareTo(this.newGeneId, o.newGeneId);
            if( r!=0 )
                return r;
            r = Utils.intsCompareTo(this.geneRgdId, o.geneRgdId);
            if( r!=0 )
                return r;
            return Utils.intsCompareTo(this.speciesTypeKey, o.speciesTypeKey);
        }
    }
}
