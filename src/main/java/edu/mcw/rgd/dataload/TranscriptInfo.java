package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: Jun 18, 2010
 * Time: 10:19:12 AM
 * represents the information extracted from a gene record about transcript and its features
 */
public class TranscriptInfo extends Transcript {

    final Log logger = LogFactory.getLog("transcript_positions");

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
