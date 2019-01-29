package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.SequenceDAO;
import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.Sequence;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Collections;
import java.util.List;

/**
 * Created by mtutaj on 11/20/2015.
 * <p>
 * for given assembly, go over all transcript acc ids and download nucleotide and protein sequences from NCBI Nucleotide database
 */
public class RefSeqSeqLoader {

    EGDAO dao = EGDAO.getInstance();

    public static void main(String[] args) throws Exception {

        boolean insertMissingSequences = true;

        new RefSeqSeqLoader().run(insertMissingSequences);
    }

    public void run(boolean insertMissingSequences) throws Exception {

        int speciesTypeKey = 3;
        int proteinSeqTypeKey = 12;

        GeneDAO geneDAO = new GeneDAO();
        SequenceDAO sequenceDAO = new SequenceDAO();
        TranscriptDAO transcriptDAO = new TranscriptDAO();

        int transcriptCount = 0;
        int nonCodingCount = 0;
        int codingCount = 0;
        int codingWithMissingProteinAccId = 0;
        int proteinSeqMissing = 0;

        List<Gene> genes = geneDAO.getAllGenes(speciesTypeKey);
        Collections.shuffle(genes);

        for (int i=0; i<genes.size(); i++ ) {
            Gene gene = genes.get(i);

            //System.out.println(i+"/"+genes.size()+". "+gene.getSymbol());

            List<Transcript> transcripts = transcriptDAO.getTranscriptsForGene(gene.getRgdId());
            for (Transcript transcript : transcripts) {
                transcriptCount++;
                if (transcript.getProteinAccId() == null) {
                    if (transcript.isNonCoding()) {
                        nonCodingCount++;
                    } else {
                        codingWithMissingProteinAccId++;
                        System.out.println("ERROR: coding transcript has no protein acc id "+ transcript.getAccId());
                    }
                    continue;
                }

                codingCount++;

                /** TODO: fix the code in the comments
                 *
                // get transcript protein as it is in NCBI database (RefSeq protein)
                Sequence seq = null;
                for (Sequence pseq : sequenceDAO.getObjectSequences(transcript.getRgdId())) {
                    if (pseq.getSeqTypeKey() != proteinSeqTypeKey) // show only RefSeq protein sequences
                        continue;
                    seq = pseq;
                    break;
                }

                if (seq == null) {
                    proteinSeqMissing++;

                    String proteinSeq = downloadProteinFasta(transcript.getProteinAccId());

                    if( insertMissingSequences ) {
                        System.out.println("INSERT protein seq for " + proteinSeq);

                        dao.createSequenceForTranscript(speciesTypeKey, proteinSeqTypeKey, transcript.getProteinAccId(),
                                transcript.getRgdId(), transcript.getProteinAccId(), proteinSeq);
                    }
                }
                 */
            }
        }

        System.out.println("transcript count  : " + transcriptCount);
        System.out.println("  non-coding count: " + nonCodingCount);
        System.out.println("  coding with missing protein acc id: " + codingWithMissingProteinAccId);
        System.out.println("  coding count    : " + codingCount);
        System.out.println("    coding with missing seq in rgd: " + proteinSeqMissing);
    }

    String downloadProteinFasta(String proteinAccId) throws Exception {

        long gi = downloadGI(proteinAccId);
        if( gi<=0 ) {
            System.out.println("Invalid GI for "+proteinAccId);
            return null;
        }

        String url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val="+gi+"&db=protein&dopt=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&maxdownloadsize=1000000";
        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(url);
        downloader.setLocalFile("data/"+proteinAccId+".fa");
        String localFile = downloader.download();

        BufferedReader reader = new BufferedReader(new FileReader(localFile));
        String line;
        StringBuilder seq = new StringBuilder();
        while( (line=reader.readLine())!=null ) {
            // skip lines starting with '>'
            if( line.startsWith(">") ) {
                continue;
            }
            seq.append(line.trim());
        }
        reader.close();

        // seq must end with '*'
        if( seq.length()>0 && seq.charAt(seq.length()-1)!='*' ) {
            seq.append("*");
        }
        return seq.toString();
    }

    long downloadGI(String proteinAccId) throws Exception {
        String url = "https://www.ncbi.nlm.nih.gov/protein/"+proteinAccId+"?report=fasta&log$=seqview&format=text";
        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(url);
        downloader.setLocalFile("data/"+proteinAccId+".xml");
        String localFile = downloader.download();

        // downloaded file is in xml format; we are looking for the following line
        // <meta name="ncbi_uidlist" content="564331920" />
        // to extract GI id
        String xml = Utils.readFileAsString(localFile);
        String pattern = "<meta name=\"ncbi_uidlist\" content=\"";
        int giPosStart = xml.indexOf(pattern);
        if( giPosStart<0 )
            return 0;
        giPosStart += pattern.length();

        int giPosEnd = xml.indexOf('\"', giPosStart);
        return Long.parseLong(xml.substring(giPosStart, giPosEnd));
    }
}
