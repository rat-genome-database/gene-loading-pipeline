package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.SequenceDAO;
import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.MappedGene;
import edu.mcw.rgd.datamodel.Sequence;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.SeqUtils;
import edu.mcw.rgd.process.Utils;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 5/27/15
 * Time: 2:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class RefSeqQcProtein {

    int genesProcessed = 0;
    int transcriptsProcessed = 0;
    int transcriptsNonCoding = 0;
    int transcriptsTripletError = 0;
    int proteinsMatchingRef = 0;
    int proteinsNotMatchingRefSameLen = 0;
    int proteinsNotMatchingRefDiffLen = 0;
    int proteinsRefSeqMissing = 0;

    public void run(int mapKey) throws Exception {

        GeneDAO geneDAO = new GeneDAO();
        SequenceDAO sequenceDAO = new SequenceDAO();
        TranscriptDAO transcriptDAO = new TranscriptDAO();


        List<MappedGene> genes = geneDAO.getActiveMappedGenes(mapKey);
        for( MappedGene gene: genes ) {
            genesProcessed++;

            List<Transcript> transcripts = transcriptDAO.getTranscriptsForGene(gene.getGene().getRgdId(), gene.getMapKey());
            for( Transcript transcript: transcripts ) {
                transcriptsProcessed++;

                if( transcript.isNonCoding() ) {
                    transcriptsNonCoding++;
                } else {
                    SeqUtils.SeqUtilsResult result = SeqUtils.translateTranscriptIntoProtein(transcript, mapKey, true);
                    if( result.tripletError ) {
                        transcriptsTripletError++;
                    } else {
                        // get transcript protein as it is in NCBI database (RefSeq protein)
                        Sequence seq = null;
                        for( Sequence pseq: sequenceDAO.getObjectSequences(transcript.getRgdId()) ) {
                            if( pseq.getSeqTypeKey()!=12 ) // show only RefSeq protein sequences
                                continue;
                            seq = pseq;
                            break;
                        }

                        if( seq==null ) {
                            proteinsRefSeqMissing++;
                            continue;
                        }

                        // translated protein always includes a stop codon '*'
                        // if a RefSeq protein does not have a stop codon, add it
                        String refSeqProtein = seq.getCloneSeq();
                        if( !refSeqProtein.endsWith("*") )
                            refSeqProtein += "*";

                        if( Utils.stringsAreEqualIgnoreCase(result.translatedProtein, refSeqProtein) ) {
                            proteinsMatchingRef++;
                        } else if( result.translatedProtein.length()==refSeqProtein.length()) {
                            proteinsNotMatchingRefSameLen++;
                        } else {
                            proteinsNotMatchingRefDiffLen++;
                        }
                    }
                }
            }
        }

        System.out.println("===== OK =====");
        System.out.println("MAP_KEY: "+mapKey);
        System.out.println("GENES PROCESSED                          : "+genesProcessed);
        System.out.println("TRANSCRIPTS PROCESSED                    : "+transcriptsProcessed);
        System.out.println("   NON-CODING                            : "+transcriptsNonCoding);
        System.out.println("   TRIPLET ERROR                         : "+transcriptsTripletError);
        System.out.println("   PROTEINS NOT MATCHING REF, SAME LENGTH: "+proteinsNotMatchingRefSameLen);
        System.out.println("   PROTEINS NOT MATCHING REF, DIFF LENGTH: "+proteinsNotMatchingRefDiffLen);
        System.out.println("   PROTEINS MATCHING REFERENCE           : "+proteinsMatchingRef);
        System.out.println("   PROTEINS REF-SEQ MISSING              : "+proteinsRefSeqMissing);
        System.out.println("===== OK =====");
    }
}
