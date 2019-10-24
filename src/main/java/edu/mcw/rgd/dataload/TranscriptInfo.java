package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * @author mtutaj
 * @since Jun 18, 2010
 * represents the information extracted from a gene record about transcript and its features
 */
public class TranscriptInfo extends Transcript {

    final Logger logger = Logger.getLogger("transcript_positions");

    private List<TranscriptLocus> loci = new ArrayList<TranscriptLocus>();

    public List<TranscriptLocus> getLoci() {
        return loci;
    }

    public TranscriptLocus createLocus(TranscriptFeature tf) throws CloneNotSupportedException {
        // look if we have locus like that one already
        for( TranscriptLocus locus: loci ) {
            if( Utils.stringsAreEqual(getAccId(), locus.getTranscriptAccId()) &&
                tf.equalsByGenomicCoords(locus.getTranscriptCoords()) ) {
                return locus;
            }
        }

        // new locus
        TranscriptLocus locus = new TranscriptLocus();
        locus.setTranscriptAccId(getAccId());
        locus.setTranscriptCoords(tf.clone());
        loci.add(locus);
        return locus;
    }

    public void updateChromosome(String chr) {
        // look if we have locus like that one already
        for( TranscriptLocus locus: loci ) {
            locus.updateChromosome(chr);
        }
    }
}
