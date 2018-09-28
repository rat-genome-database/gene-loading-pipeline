package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.Sequence;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;

import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 10/29/12
 * Time: 8:52 AM
 */
public class RefSeqLoader {

    private List<String> ncbiFiles;

    int countConflictMultiTranscript = 0;
    int countConflictNoTranscriptMatch = 0;
    int countSequenceMatch = 0;
    int countSequenceUpdate = 0;
    int countSequenceNew = 0;

    public void run() throws Exception {

        // go through every file
        int speciesTypeKey;
        boolean isProtein = false; // sequence is either protein or rna

        // go through every file
        for(String fileName: ncbiFiles) {

            // detect species from file name
            if( fileName.contains("rat") )
                speciesTypeKey = SpeciesType.RAT;
            else if( fileName.contains("mouse") )
                speciesTypeKey = SpeciesType.MOUSE;
            else if( fileName.contains("human") )
                speciesTypeKey = SpeciesType.HUMAN;
            else {
                System.out.println("Unknown species for file named "+fileName);
                continue;
            }

            // detect sequence type (protein or rna) from file name
            if( fileName.contains("protein") )
                isProtein = true;
            else if( fileName.contains("rna") )
                isProtein = false;
            else {
                System.out.println("Unknown sequence type (neither protein or rna) for file named "+fileName);
                continue;
            }

            // get the file
            FileDownloader fd = new FileDownloader();
            fd.setExternalFile(fileName);
            fd.setLocalFile("data/"+SpeciesType.getCommonName(speciesTypeKey).toLowerCase()+(isProtein?"_protein":"_rna")+".faa.gz");
            System.out.println("downloading file "+fileName);
            fd.setPrependDateStamp(true);
            String localFile = fd.downloadNew();
            System.out.println("downloaded "+localFile);

            //String localFile = isProtein ? "/tmp/refseq/rat.protein.faa.gz" : "/tmp/refseq/rat.rna.fna.gz";
            BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(localFile))));

            String line, gi = "", accId = "", info = "", seqData = "";
            while( (line = reader.readLine())!= null ) {

                if( line.startsWith(">") ) {

                    // save sequence to database
                    save(isProtein, speciesTypeKey, gi, accId, info, seqData);
                    gi = accId = info = seqData = "";

                    String[] cols = line.split("[|]", -1);
                    if( !cols[0].equals(">gi") || !cols[2].equals("ref") ) {
                        System.out.println("bad input line for a transcript: "+line);
                        continue;
                    }
                    gi = cols[1];
                    accId = cols[3];
                    info = cols[4];
                }
                else {
                    seqData += line.trim();
                }
            }

            reader.close();
        }

        // display counters
        System.out.println("conflict accession ids assigned to multiple transcripts: "+countConflictMultiTranscript);
        System.out.println("conflict accession id does not match a transcript: "+countConflictNoTranscriptMatch);
        System.out.println("sequences matched in RGD: "+countSequenceMatch);
        System.out.println("sequences updated in RGD: "+countSequenceUpdate);
        System.out.println("sequences inserted into RGD: "+countSequenceNew);
    }

    void save(boolean isProtein, int speciesTypeKey, String gi, String accId, String info, String data) throws Exception {

        if( accId==null || accId.isEmpty() )
            return;
        EGDAO dao = EGDAO.getInstance();

        // strip version part of accession id (XM_216369.3' ==> 'XM_216369')
        int dotPos = accId.lastIndexOf('.');
        if( dotPos > 0 )
            accId = accId.substring(0, dotPos);

        // convert transcript accession id into transcript rgd id
        List<Transcript> transcripts = isProtein ? dao.getTranscriptsByProteinAccId(accId) : dao.getTranscriptsByAccId(accId);

        validateAgainstActiveGenes(transcripts, dao);

        if( transcripts.size()>1 ) {
            System.out.println("CONFLICT: acc_id "+accId+" is assigned to multiple transcripts");
            countConflictMultiTranscript++;
            return;
        }
        if( transcripts.isEmpty() ) {
            System.out.println("CONFLICT: acc_id "+accId+" is not assigned to a transcript");
            countConflictNoTranscriptMatch++;
            return;
        }
        Transcript transcript = transcripts.get(0);

        // get existing sequences for the transcript
        List<Sequence> seqsInRgd = dao.getObjectSequences(transcript.getRgdId());
        int seqTypeKey = isProtein ? 12 : 7;

        // look for sequence of type mRna or protein
        for( Sequence seq: seqsInRgd ) {
            if( seq.getSeqTypeKey()==seqTypeKey ) {
                // mRna sequence has been found -- compare the seq data
                if( Utils.stringsAreEqualIgnoreCase(data, seq.getCloneSeq()) ) {
                    // sequences are the same
                    dao.updateLastModifiedDate(seq.getRgdId());
                    countSequenceMatch++;
                }
                else {
                    dao.updateSeqData(seq.getSeqKey(), data);
                    dao.updateLastModifiedDate(seq.getRgdId());
                    //if( isProtein )
                    //    System.out.println("SEQ-UPDATE "+transcript.getProteinAccId());
                    //else
                    //    System.out.println("SEQ-UPDATE "+transcript.getAccId());
                    countSequenceUpdate++;
                }
                return; // sequence has been processed!
            }
        }

        // sequence not found -- create a new sequence
        dao.createSequenceForTranscript(speciesTypeKey, seqTypeKey, info, transcript.getRgdId(), accId, data);
        //if( isProtein )
        //  System.out.println("SEQ-INSERT "+transcript.getProteinAccId());
        //else
        //    System.out.println("SEQ-INSERT "+transcript.getAccId());
        countSequenceNew++;
    }

    void validateAgainstActiveGenes(List<Transcript> transcripts, EGDAO dao) throws Exception {

        Iterator<Transcript> it = transcripts.iterator();
        while( it.hasNext() ) {
            Transcript tr = it.next();
            RgdId rgdId = dao.getRgdId(tr.getGeneRgdId());
            if( rgdId==null ) {
                it.remove();
            }
            if( !Utils.stringsAreEqual(rgdId.getObjectStatus(),"ACTIVE") ) {
                it.remove();
            }
        }
    }

    public void setNcbiFiles(List<String> ncbiFiles) {
        this.ncbiFiles = ncbiFiles;
    }

    public List<String> getNcbiFiles() {
        return ncbiFiles;
    }

}
