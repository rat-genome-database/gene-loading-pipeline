package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.FileDownloader;

import java.io.*;
import java.net.*;
import java.sql.ResultSet;
import java.util.*;
import java.util.zip.*;


/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: Mar 29, 2010
 * Time: 5:20:09 PM
 * used to check all RefSeq references from NCBI
 */
public class RefSeqValidator {

    // list of all input files to be checked
    List<String> ncbiFiles;

    // ftp directory with files containing removed RefSeq records (ACC_IDs that have been replaced or removed)
    String removedRecordsDir;

    // summaries
    int refSeqTotal = 0; // total number of RefSeq records analyzed
    int rgdIdsTotal = 0; // total number of ACC_ID-RGD_ID pairs analyzed
    int refSeqMissing = 0; // not found in NCBI RefSeq, found in Rgd
    int rgdIdNotActive = 0;
    int diffSpecies = 0;
    int wrongXdbKey = 0;
    int rgdMissing = 0;
    int removedRecordsTotal = 0; // total number of records read from /removed/ folder
    int mismatchesAnnotated = 0; // number of mismatched lines that have been annotated from removed-records

    public void run() throws Exception {

        // maps the locus n
        int initialCapacity = 360007; // we handle over 300k records so it makes sense to preset hashmap to decent capacity 
        Map<String,RefSeqInfo> map = new HashMap<String,RefSeqInfo>(initialCapacity);

        // go through every file
        for(String fileName: ncbiFiles) {

            // get the file through ftp
            URL url = new URL(fileName+";type=i");
            URLConnection urlc = url.openConnection();
            InputStream is = urlc.getInputStream(); // To download

            // the file is supposed to be a character file
            BufferedReader reader;
            if(fileName.endsWith("gz")) {
                reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(is)));
            }
            else {
                reader = new BufferedReader(new InputStreamReader(is));
            }

            // determine species and whether it is a protein accession
            int species = fileName.contains("norvegicus") ? SpeciesType.RAT :
                    fileName.contains("musculus") ? SpeciesType.MOUSE :
                    fileName.contains("sapiens") ? SpeciesType.HUMAN :
                    SpeciesType.ALL;
            boolean isProtein = fileName.contains("gpff");
            RefSeqInfo rs;

            String line, accid;
            while( (line = reader.readLine())!= null ) {
                if( line.startsWith("ACCESSION")) {
                    // 2nd and other words are ACCESSION ACC_ID
                    String[] words = line.split("\\s+");
                    if( words!=null && words.length>1 ) {
                        for( int i=1; i<words.length; i++ ) {
                            accid = words[i];
                            rs = map.get(accid);
                            if( rs == null ) {
                                rs = new RefSeqInfo(accid, species, isProtein);
                                map.put(accid, rs);
                            }
                            //System.out.println(rs);
                        }
                    }
                }
            }

            reader.close();
        }

        this.refSeqTotal = map.size();
        System.out.println("ref-seq count: "+refSeqTotal);

        // get all accession_ids from database, and update our map
        matchRefSeqs(map);

        // process removed RefSeq records
        processRemovedRecords(map);

        // print mismatched records
        printMismatches(map);

        // print stats
        printStats();
    }

    // match all records from these in database
    void matchRefSeqs(Map<String,RefSeqInfo> map) throws Exception {

        ResultSet rs = EGDAO.getInstance().getAllXdbIds();
        while( rs.next() ) {
            // access fields from result set
            String accId = rs.getString(1);
            int xdbKey = rs.getInt(2);
            int rgdId = rs.getInt(3);
            String status = rs.getString(4);
            int species = rs.getInt(5);

            rgdIdsTotal++;
            //System.out.printf("%d. %s RGD_ID=%d\n", rgdIdsTotal, accId, rgdId);

            // do the checks
            RefSeqInfo info = map.get(accId);
            if( info==null ) {
                // not found in NCBI
                info = new RefSeqInfo(accId, species, false);
                info.addRgdIdEntry(rgdId, status, xdbKey, "not found in NCBI, found in RGD, ");
                map.put(accId, info);
                this.refSeqMissing++;
                continue;
            }

            // found in both NCBI and RGD
            String error = "";
            if( !status.equals("ACTIVE") ) {
                error += "RGD_ID is not active, ";
                rgdIdNotActive++;
            }
            if( species!=info.species ) {
                error += "different species in RGD, ";
                diffSpecies++;
            }
            if( xdbKey!=1 && xdbKey!=7 && xdbKey!=10 ) {
                error += "XDB_KEY not (1,7,10), ";
                wrongXdbKey++;
            }
            info.addRgdIdEntry(rgdId, status, xdbKey, error);
        }

        // post processing: mark as error records having empty RGD_ID
        for( Map.Entry entry: map.entrySet() ) {
            RefSeqInfo info = (RefSeqInfo) entry.getValue();
            if( info.rgdIdEntries==null ) {
                info.error += "not found in RGD, ";
                rgdMissing ++;
            }
        }
        System.out.println("total count of ref-seq records: "+map.size());

        // remove from the map all records with perfect matches
        Iterator it = map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry entry = (Map.Entry)it.next();
            RefSeqInfo info = (RefSeqInfo) entry.getValue();
            // remove all rgd records with empty error message
            if( info.rgdIdEntries!=null ) {
                Iterator it2 = info.rgdIdEntries.iterator();
                while(it2.hasNext() ) {
                    RefSeqInfo.RgdIdEntry rentry = (RefSeqInfo.RgdIdEntry) it2.next();
                    if( rentry.error==null || rentry.error.length()==0 ) {
                        it2.remove();
                    }
                }
                if( info.rgdIdEntries.size()==0 )
                    info.rgdIdEntries = null;
            }

            // remove all map entries having no mismatched rgdids and no error flag set
            if( info.rgdIdEntries==null && (info.error==null || info.error.length()==0) ) {
                it.remove();
            }
        }

        System.out.println("total count of mismatched ref-seq records: "+map.size());
    }

    void printMismatches(Map<String,RefSeqInfo> map) throws Exception {

        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("data/refseq_mismatches.tcv")));
        out.format("Species\tError-message\tRGD-ID\tACC-ID\tStatus\tXDB-KEY\tNotes\n");
        for( Map.Entry entry: map.entrySet() ) {
            RefSeqInfo info = (RefSeqInfo) entry.getValue();
            String annotation = info.getRemovedRecords();
            if( info.error.length()>0 ) {
                out.format("%d\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    info.species, info.error, "", info.accId, "", "", annotation );

                if( annotation.length()>0 )
                    this.mismatchesAnnotated++;
            }
            if( info.rgdIdEntries==null )
                continue;
            for( RefSeqInfo.RgdIdEntry rentry: info.rgdIdEntries ) {
                if( rentry.error!=null && rentry.error.length()>0 ) {
                    out.format("%d\t%s\t%d\t%s\t%s\t%d\t%s\n",
                        info.species, info.error+rentry.error, rentry.rgdId, info.accId,
                        rentry.status, rentry.xdbKey, annotation );

                    if( annotation.length()>0 )
                        this.mismatchesAnnotated++;
                }
            }
        }
        out.close();
    }

    void printStats() throws Exception {

        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("data/refseq_stats.tcv")));

        out.println("total number of RefSeq records analyzed: "+refSeqTotal);
        out.println("total number of ACC_ID-RGD_ID pairs analyzed: "+rgdIdsTotal);
        out.println("not found in NCBI RefSeq, but found in RGD: "+refSeqMissing);
        out.println("RGD_ID not active: "+rgdIdNotActive);
        out.println("species mismatch: "+diffSpecies);
        out.println("XDB_KEY not 1,7, or 10: "+wrongXdbKey);
        out.println("ACC_ID not found in RGD: "+rgdMissing);
        out.println("Number of lines read from "+removedRecordsDir+": "+removedRecordsTotal);
        out.println("Number of mismatches annotated through removed-records: "+mismatchesAnnotated);
        out.close();
    }

    // this methods loads all removed refseqs and sorts them by accId and date
    void processRemovedRecords(Map<String,RefSeqInfo> mapMismatches) throws Exception {

        // go to source directory and do FTP listing
        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(removedRecordsDir);
        String[] files = downloader.listFiles();

        // sort the loaded files by date
        Arrays.sort(files, new Comparator<String>() {
            public int compare(String f1, String f2) {
                //                                            0123456789012345678901234
                // the file name should have at the end date: removed-records.mmdd.yyyy
                // change it to:                                              yyyy.mmdd
                // we want files to be sorted by date!
                String fname1 = f1==null ? "" : f1;
                String fname2 = f2==null ? "" : f2;
                if( fname1.length()!=25 || fname2.length()!=25 )
                    return fname1.compareTo(fname2); // file name to short

                String pattern1 = fname1.substring(21,25)+fname1.substring(16,20);
                String pattern2 = fname2.substring(21,25)+fname2.substring(16,20);
                return pattern1.compareTo(pattern2);
            }
        });

        for( String file: files ) {
            downloader.setExternalFile(removedRecordsDir+'/'+file);
            //System.out.println("loading "+file);
            downloader.setLocalFile("data/"+file);

            // optimization: verify if the local file already exists to avoid downloading
            String localFile = downloader.getLocalFile();
            File rFile = new File(localFile);
            if( !rFile.exists() )
                localFile = downloader.download();

            BufferedReader reader = new BufferedReader(new FileReader(localFile));
            String line, accId;
            while( (line=reader.readLine())!=null ) {
                removedRecordsTotal++;

                // columns are separated by tab chars;
                // col 0: ACC_ID, col 3: message, col 4: date
                String[] words = line.split("\\t");
                if( words.length>=5 ) {
                    accId = words[0];
                    RefSeqInfo ri = mapMismatches.get(accId);
                    if( ri!=null ) {
                        ri.addRemovedRecordMsg(words[4], words[3]);
                    }
                }
            }
            reader.close();
        }
    }

    public List<String> getNcbiFiles() {
        return ncbiFiles;
    }

    public void setNcbiFiles(List<String> ncbiFiles) {
        this.ncbiFiles = ncbiFiles;
    }

    public String getRemovedRecordsDir() {
        return removedRecordsDir;
    }

    public void setRemovedRecordsDir(String removedRecordsDir) {
        this.removedRecordsDir = removedRecordsDir;
    }


    // information for a particular refseq
    class RefSeqInfo {
        public String accId;
        public int species; // species type key
        public boolean protein; // protein (gpff file) or mrna (gbff)
        public String error;
        public ArrayList<RgdIdEntry> rgdIdEntries;
        public String removedRecordMsg; // list of messages about removed refseq records

        public void addRemovedRecordMsg(String date, String msg) {
            if( removedRecordMsg==null )
                removedRecordMsg = '\t' + date+": "+msg;
            else
                removedRecordMsg += '\t' + date+": "+msg;
        }

        // return textual representation of removed records
        public String getRemovedRecords() {
            if( removedRecordMsg==null )
                return "";
            return removedRecordMsg;
        }

        public RefSeqInfo(String accId, int species, boolean protein) {
            this.accId = accId;
            this.species = species;
            this.protein = protein;
            this.error = "";
        }


        @Override
        public String toString() {
            return String.format("%s - %s", accId, SpeciesType.getCommonName(species));
        }

        public RgdIdEntry addRgdIdEntry(int rgdId, String status, int xdbKey, String error) {
            if(rgdIdEntries==null) {
                rgdIdEntries = new ArrayList<RgdIdEntry>();
            }
            RgdIdEntry entry = new RgdIdEntry();
            entry.rgdId = rgdId;
            entry.xdbKey = xdbKey;
            entry.status = status;
            entry.error = error;
            rgdIdEntries.add(entry);
            return entry;
        }

        public class RgdIdEntry {
            public int rgdId;   // RGD_ID associated with this accession field; 0 means this record is not found in RGD_ACC_XDB table
            public int xdbKey; // XDB_KEY accociated with this ACC_ID
            public String status; // as retrieved from RGD_IDS table; should be ACTIVE
            public String error;
        }
    }
}
