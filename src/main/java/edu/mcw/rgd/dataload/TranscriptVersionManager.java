package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.process.CounterPool;

import java.util.HashMap;
import java.util.Map;

/**
 * manages transcript versions in STABLE_TRANSCRIPTS table
 */
public class TranscriptVersionManager {

    // singleton
    private static TranscriptVersionManager _instance = new TranscriptVersionManager();

    public static TranscriptVersionManager getInstance() {
        return _instance;
    }

    private TranscriptVersionManager() {
    }

    // hashmap of tr accession to Info
    private Map<String, Info> map = new HashMap<String, Info>();

    public void addVersion(String acc, String version) {

        int versionNr = Integer.parseInt(version);

        Info info = map.get(acc);
        if( info==null ) {
            info = new Info();
            info.version = versionNr;
            map.put(acc, info);
        } else {
            if( info.version==0 ) {
                info.version = versionNr;
            }
            else if( info.version != versionNr ) {
                throw new RuntimeException("transcript acc ver mismatch: "+info.version+" vs "+versionNr);
            }
        }
    }

    public void addRgdId(String acc, int rgdId) {

        Info info = map.get(acc);
        if( info==null ) {
            info = new Info();
            info.rgdId = rgdId;
            map.put(acc, info);
        } else {
            if( info.rgdId==0 ) {
                info.rgdId = rgdId;
            }
            else if( info.rgdId != rgdId ) {
                throw new RuntimeException("transcript rgd id mismatch: "+info.rgdId+" vs "+rgdId);
            }
        }
    }

    public void qcAndLoad(CounterPool counters) throws Exception {

        EGDAO dao = EGDAO.getInstance();

        for( Map.Entry<String, Info> entry: map.entrySet() ) {
            String acc = entry.getKey();
            Info info = entry.getValue();
            counters.increment("TRANSCRIPT_VERSIONS__PROCESSED");
            String accVer = acc + "." + info.version;

            String accVerInDb = dao.getTranscriptVersionInfo(acc);
            if( accVerInDb==null ) {
                dao.insertTranscriptVersionInfo(acc, accVer, info.rgdId);
                counters.increment("TRANSCRIPT_VERSIONS_INSERTED");
            } else {
                if( accVerInDb.equals(accVer) ) {
                    counters.increment("TRANSCRIPT_VERSIONS_UP_TO_DATE");
                } else {
                    dao.updateTranscriptVersionInfo(acc, accVer);
                    counters.increment("TRANSCRIPT_VERSIONS_MODIFIED");
                }
            }
        }
    }

    class Info {
        public int version;
        public int rgdId;

        public String toString() {
            return "VER="+version+", RGD:"+rgdId;
        }
    }
}
