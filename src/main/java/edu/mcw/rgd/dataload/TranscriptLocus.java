package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.*;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * @author mtutaj
 * @since Mar 19, 2012
 * represents the information extracted from a gene record about a transcript locus and its features;
 * note: one transcript could have multiple loci on one assembly map
 */
public class TranscriptLocus {

    String transcriptAccId; // if known
    TranscriptFeature transcriptCoords; // coordinates for transcript
    TranscriptFeature utr5; // coordinates of 5 prime UTR, if any
    TranscriptFeature utr3; // coordinates of 3 prime UTR, if any

    List<TranscriptFeature> coords; // genomic coordinates -- exons
    List<TranscriptFeature> cdsCoords; // genomic coordinates -- cdss
    // note: we will cheat a bit here, and put into bandType field the type of genomic feature: Feature


    public TranscriptLocus() {
        coords = new ArrayList<TranscriptFeature>();
        cdsCoords = new ArrayList<TranscriptFeature>();
    }

    public void addGenomicCoords(TranscriptFeature md) {
        coords.add(md);
    }

    // adjust transcript pos based on first and last exon
    public void adjustTranscriptPos() {
        int startPos = -1, stopPos = -1;
        for( TranscriptFeature exon: coords ) {
            // new startPos?
            if( startPos<0 || exon.getStartPos()<startPos ) {
                startPos = exon.getStartPos();
            }
            // new stopPos?
            if( stopPos<0 || exon.getStopPos()>stopPos ) {
                stopPos = exon.getStopPos();
            }
        }

        // do we have a new position available?
        if( startPos<0 || stopPos<0 ) {
            return; // no position available based on exons
        }
        // see if transcript pos changed
        if( transcriptCoords.getStartPos()!=startPos || transcriptCoords.getStopPos()!=stopPos ) {

            final Logger log = Logger.getLogger("process");
            log.debug("Transcript "+transcriptAccId+" pos adjusted: "+
                "["+transcriptCoords.getStartPos()+", "+transcriptCoords.getStopPos()+"] ==> "+
                "["+startPos+", "+stopPos+"]");

            transcriptCoords.setStartPos(startPos);
            transcriptCoords.setStopPos(stopPos);
        }
    }

    public void addGenomicCoordsForPeptide(TranscriptFeature md) {
        cdsCoords.add(md);
    }

    public void buildUtrs() {

        // check if the transcript is non-coding
        if( this.cdsCoords.size()==0 ) {
            //this.setNonCoding(true);
            return; // there will be no UTR's for non-coding regions
        }

        // is there anything to do?
        if( this.coords.size()==0 )
            return;
        String strand = this.coords.get(0).getStrand();
        String chromosome = this.coords.get(0).getChromosome();
        int mapKey = this.coords.get(0).getMapKey();

        // build arrays of exon coordinates and cds coordinates
        List<CoordPos> exonPositions = new ArrayList<CoordPos>(this.coords.size()*2);
        List<CoordPos> cdsPositions = new ArrayList<CoordPos>(this.cdsCoords.size()*2);
        for( MapData md: this.coords ) {
            exonPositions.add(new CoordPos(md.getStartPos(), true));
            exonPositions.add(new CoordPos(md.getStopPos(), false));
        }
        for( MapData md: this.cdsCoords ) {
            cdsPositions.add(new CoordPos(md.getStartPos(), true));
            cdsPositions.add(new CoordPos(md.getStopPos(), false));
        }
        Collections.sort(exonPositions);
        Collections.sort(cdsPositions);

        // determine stop of 1st UTR region,
        Iterator<CoordPos> itExon = exonPositions.iterator();
        CoordPos exonCoord = itExon.next();
        int pos = exonCoord.getPos();
        if( cdsPositions.size()==0 ) {
            // there is no CDS region at all -- UTR only
            TranscriptFeature md = new TranscriptFeature();
            md.setStartPos(pos);
            md.setStopPos(exonPositions.get(exonPositions.size()-1).getPos());
            md.setChromosome(chromosome);
            md.setStrand(strand);
            md.setMapKey(mapKey);
            addFirstUtr(md);
            return;
        }

        // there is a CDS region - does it have a first utr
        int cdsPos = cdsPositions.get(0).getPos();
        if( cdsPos > pos ) {
            // there is first utr region!
            TranscriptFeature md = new TranscriptFeature();
            md.setStartPos(pos);
            md.setStopPos(cdsPos-1);
            md.setChromosome(chromosome);
            md.setStrand(strand);
            md.setMapKey(mapKey);
            addFirstUtr(md);
        }

        // read last exon and cds stop positions
        pos = exonPositions.get(exonPositions.size()-1).getPos();
        cdsPos = cdsPositions.get(cdsPositions.size()-1).getPos();
        if( cdsPos < pos ) {
            // there is last utr region!
            TranscriptFeature md = new TranscriptFeature();
            md.setStartPos(cdsPos+1);
            md.setStopPos(pos);
            md.setChromosome(chromosome);
            md.setStrand(strand);
            md.setMapKey(mapKey);
            addLastUtr(md);
        }
    }

    void addFirstUtr(TranscriptFeature md) {
        if( md.getStrand().equals("+") ) {
            // positive strand starts from 5'UTR
            md.setFeatureType(TranscriptFeature.FeatureType.UTR5);
            this.utr5 = md;
        }
        else {
            // negative strand starts with 3'UTR
            md.setFeatureType(TranscriptFeature.FeatureType.UTR3);
            this.utr3 = md;
        }
    }

    void addLastUtr(TranscriptFeature md) {
        if( md.getStrand().equals("+") ) {
            // positive strand ends with 3'UTR
            md.setFeatureType(TranscriptFeature.FeatureType.UTR3);
            this.utr3 = md;
        }
        else {
            // negative strand ends with 5'UTR
            md.setFeatureType(TranscriptFeature.FeatureType.UTR5);
            this.utr5 = md;
        }
    }

    public List<TranscriptFeature> getCoords() {
        return coords;
    }

    public void setCoords(List<TranscriptFeature> coords) {
        this.coords = coords;
    }

    public TranscriptFeature getTranscriptCoords() {
        return transcriptCoords;
    }

    public void setTranscriptCoords(TranscriptFeature transcriptCoords) {
        this.transcriptCoords = transcriptCoords;
    }

    public String getTranscriptAccId() {
        return transcriptAccId;
    }

    public void setTranscriptAccId(String transcriptAccId) {
        this.transcriptAccId = transcriptAccId;
    }

    class CoordPos implements Comparable<CoordPos> {
        int pos;
        public boolean isStart;

        public CoordPos(int pos, boolean isStart) {
            this.pos = pos;
            this.isStart = isStart;
        }

        public int compareTo(CoordPos o) {
            return this.pos - o.pos;
        }

        public int getPos() {
            return this.pos;
        }
    }

    public TranscriptFeature getUtr5() {
        return utr5;
    }

    public void setUtr5(TranscriptFeature utr5) {
        this.utr5 = utr5;
    }

    public TranscriptFeature getUtr3() {
        return utr3;
    }

    public void setUtr3(TranscriptFeature utr3) {
        this.utr3 = utr3;
    }

    public void updateChromosome(String chromosome) {

        if( transcriptCoords.getChromosome()==null )
            transcriptCoords.setChromosome(chromosome);
        else
            chromosome = transcriptCoords.getChromosome();

        for( TranscriptFeature tf: coords ) {
            tf.setChromosome(chromosome);
        }
        for( TranscriptFeature tf: cdsCoords ) {
            tf.setChromosome(chromosome);
        }

        this.buildUtrs();

        if( utr5!=null ) {
            utr5.setChromosome(chromosome);
        }
        if( utr3!=null ) {
            utr3.setChromosome(chromosome);
        }
    }
}
