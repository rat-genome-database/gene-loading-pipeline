package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

public class LoadTranscriptsFromGff3 {

    public static void main(String[] args) throws IOException {

        try {
            new LoadTranscriptsFromGff3().run();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    final int mapKey = 1311;
    BulkGeneLoaderImpl impl = new BulkGeneLoaderImpl();

    Map<String, String> chrMap = new HashMap<>();

    void run() throws Exception {

        String fname = "/Users/mtutaj/Downloads/r/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.gff";

        MapDAO mapDAO = new MapDAO();
        List<Chromosome> chrList = mapDAO.getChromosomes(mapKey);
        for( Chromosome chr: chrList ) {
            chrMap.put(chr.getRefseqId(), chr.getChromosome());
        }

        // gene map: ncbi-gene-id -> list of TrInfo
        Map<String, List<TrInfo>> geneMap = loadGeneMap(fname);

        for( Map.Entry<String, List<TrInfo>> entry: geneMap.entrySet() ) {
            processGene(entry.getKey(), entry.getValue());
        }
    }

    void processGene(String egId, List<TrInfo> transcripts ) throws Exception {
        // match
        BulkGene bg = loadBulkGene(egId, transcripts);

        impl.handleTranscripts(bg);
    }

    BulkGene loadBulkGene(String egId, List<TrInfo> transcripts) throws Exception {

        BulkGene bg = new BulkGene();
        bg.setEgId(egId);

        // resolve bg id
        EGDAO dao = EGDAO.getInstance();
        List<Gene> genes = dao.getGenesByEGID(egId);
        if( genes.isEmpty() ) {
            System.out.println("NO GENES");
        } else if( genes.size()>1 ) {
            System.out.println("MULTI GENES");
        } else {
            int geneRgdId = genes.get(0).getRgdId();

            Flags flags = new Flags();
            flags.setRgdId(geneRgdId);
            bg.setCustomFlags(flags);

            bg.rgdTranscripts = dao.getNcbiTranscriptsForGene(geneRgdId);

            // incoming transcripts
            bg.transcripts = new LinkedList<>();
            for( TrInfo trInfo: transcripts ) {
                TranscriptInfo t = new TranscriptInfo();
                t.setGeneRgdId(geneRgdId);
                // strip tr version
                String acc = trInfo.acc;
                int dotPos = acc.lastIndexOf('.');
                String accId = dotPos>0 ? acc.substring(0, dotPos) : acc;
                t.setAccId(accId);
                MapData md = new MapData();
                md.setMapKey(mapKey);
                md.setStrand(trInfo.strand);
                md.setStartPos(trInfo.startPos);
                md.setStopPos(trInfo.stopPos);
                md.setChromosome(chrMap.get(trInfo.chrAcc));
                t.getGenomicPositions().add(md);
                bg.transcripts.add(t);
            }
        }

        return bg;
    }

    Map<String, List<TrInfo>> loadGeneMap(String fname) throws IOException {

        BufferedReader in = Utils.openReader(fname);
        String line;
        Map<String, Integer> objCount = new HashMap<>();
        Map<String, Integer> objCountSkipped = new HashMap<>();

        Map<String, List<TrInfo>> geneMap = new HashMap<>();

        while( (line=in.readLine())!=null ) {
            if( geneMap.size()==10 ) break;

            // skip comment lines
            if( line.startsWith("#") ) {
                continue;
            }
            // valid gff line must have 9 columns
            String[] cols = line.split("[\\t]", -1);
            if( cols.length<9 ) {
                continue;
            }

            String chrAcc = cols[0]; // f.e.NC_023642.1
            String db = cols[1]; // f.e. Gnomon
            String obj = cols[2]; // f.e. gene
            int startPos = Integer.parseInt(cols[3]);
            int stopPos = Integer.parseInt(cols[4]);
            String dot = cols[5];
            String strand = cols[6];
            String phase = cols[7];
            String info = cols[8];

            Integer cnt = objCount.get(obj);
            if (cnt == null) cnt = 1;
            else cnt++;
            objCount.put(obj, cnt);

            if( obj.equals("mRNA") || obj.equals("lnc_RNA") || obj.equals("transcript") || obj.equals("tRNA") || obj.equals("rRNA") ) {
                String geneId = getTokenValue(info, "Dbxref=GeneID:", ",");
                String trAcc = getTokenValue(info, "Genbank:", ";");
                String trId = getTokenValue(info, "ID=", ";");

                if( geneId!=null && trAcc!=null ) {
                    List<TrInfo> list = geneMap.get(geneId);
                    if( list == null ) {
                        list = new ArrayList<>();
                        geneMap.put(geneId, list);
                    }
                    TrInfo trInfo = null;
                    for( TrInfo i: list ) {
                        if( i.id.equals(trId) ) {
                            trInfo = i;
                            break;
                        }
                    }
                    if( trInfo == null ) {
                        trInfo = new TrInfo();
                        trInfo.id = trId;
                        trInfo.acc = trAcc;
                        trInfo.chrAcc = chrAcc;
                        trInfo.strand = strand;
                        trInfo.startPos = startPos;
                        trInfo.stopPos = stopPos;
                        list.add(trInfo);
                    }
                }
                continue;
            }

            if( obj.equals("exon") ) {
                String geneId = getTokenValue(info, "Dbxref=GeneID:", ",");
                String trId = getTokenValue(info, "Parent=", ";");

                if( geneId!=null && trId!=null ) {
                    List<TrInfo> list = geneMap.get(geneId);
                    if( list != null ) {
                        TrInfo trInfo = null;
                        for( TrInfo i: list ) {
                            if( i.id.equals(trId) ) {
                                trInfo = i;
                                break;
                            }
                        }

                        ExonInfo e = new ExonInfo();
                        e.startPos = startPos;
                        e.stopPos = stopPos;
                        e.strand = strand;
                        trInfo.exons.add(e);
                    }
                }
                continue;
            }

            if( obj.equals("CDS") ) {
                String geneId = getTokenValue(info, "Dbxref=GeneID:", ",");
                String trId = getTokenValue(info, "Parent=", ";");

                if( geneId!=null && trId!=null ) {
                    List<TrInfo> list = geneMap.get(geneId);
                    if( list == null ) {
                        System.out.println("unexpected 2");
                    }
                    TrInfo trInfo = null;
                    for( TrInfo i: list ) {
                        if( i.id.equals(trId) ) {
                            trInfo = i;
                            break;
                        }
                    }

                    if( trInfo.cdsStart==0 ) {
                        trInfo.cdsStart = startPos;
                    } else if( startPos < trInfo.cdsStart ) {
                        trInfo.cdsStart = startPos;
                    }

                    if( trInfo.cdsStop==0 ) {
                        trInfo.cdsStop = stopPos;
                    } else if( stopPos > trInfo.cdsStop ) {
                        trInfo.cdsStop = stopPos;
                    }
                }
                continue;
            }

            cnt = objCountSkipped.get(obj);
            if (cnt == null) cnt = 1;
            else cnt++;
            objCountSkipped.put(obj, cnt);
        }
        in.close();

        System.out.println("objCount: ");
        for( Map.Entry<String, Integer> entry: objCount.entrySet() ) {
            System.out.println("    "+entry.getKey()+": "+entry.getValue());
        }

        System.out.println("objCountSkipped: ");
        for( Map.Entry<String, Integer> entry: objCountSkipped.entrySet() ) {
            System.out.println("    "+entry.getKey()+": "+entry.getValue());
        }
        return geneMap;
    }

    String getTokenValue(String info, String startToken, String endToken) {
        String result = null;
        int p1 = info.indexOf(startToken);
        if( p1>=0 ) {
            int p2 = info.indexOf(endToken, p1);
            if( p2>=0 ) {
                result = info.substring(p1 + startToken.length(), p2);
            }
        }
        return result;
    }

    static public class TrInfo {
        String id;  // rna-XM_008004389.1, etc
        String acc; // XM_xxx etc
        String chrAcc;
        String strand; // '+' or '-'
        int startPos;
        int stopPos;

        int cdsStart;
        int cdsStop;
        List<ExonInfo> exons = new ArrayList<>();
    }

    static public class ExonInfo {
        String strand; // '+' or '-'
        int startPos;
        int stopPos;
    }
}
