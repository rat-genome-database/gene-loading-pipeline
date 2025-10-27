package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.CounterPool;
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

    final int mapKey = 306;
    CounterPool counters = new CounterPool();

    Map<String, String> chrMap = new HashMap<>();

    void run() throws Exception {

        String fname = "/Users/mtutaj/Downloads/r/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.gff";
        fname = "/tmp/g/306.gff";

        MapDAO mapDAO = new MapDAO();
        List<Chromosome> chrList = mapDAO.getChromosomes(mapKey);
        for( Chromosome chr: chrList ) {
            chrMap.put(chr.getRefseqId(), chr.getChromosome());
        }

        // gene map: ncbi-gene-id -> list of GeneInfo
        Map<String, GeneInfo> geneMap = loadGeneMap(fname);

        List<GeneInfo> randomizedList = new ArrayList<>(geneMap.values());
        Collections.shuffle(randomizedList);

        for( GeneInfo geneInfo: randomizedList ) {
            processGene(geneInfo);
        }

        System.out.println(counters.dumpAlphabetically());
    }

    Map<String, GeneInfo> loadGeneMap(String fname) throws Exception {

        BufferedReader in = Utils.openReader(fname);
        String line;
        Map<String, Integer> objCount = new HashMap<>();

        Map<String, GeneInfo> geneMap = new HashMap<>();
        String chr = "", regionChrAcc = "";
        HashSet<String> ignoredFeatures = new HashSet<>(); // we skip these features together with their exon child objects

        int lineNr = 0;

        while( (line=in.readLine())!=null ) {

            lineNr++;

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

            switch (obj) {
                case "region" -> {
                    regionChrAcc = chrAcc;
                    chr = getTokenValue(info, "Name=", ";");
                    if( chr!=null ) {
                        System.out.println("processing chromosome " + chr);
                    }
                }

                case "pseudogene", "gene" -> {
                    String ncbiGeneId = getTokenValue(info, "Dbxref=GeneID:", ",", ";");
                    String geneSymbol = getTokenValue(info, "Name=", ";");
                    String geneBioType = getTokenValue(info, "gene_biotype=", ";");
                    String pseudoStr = getTokenValue(info, "pseudo=");
                    String geneRgdIdStr = getTokenValue(info, ",RGD:", ";");
                    int geneRgdId = 0;
                    if( !Utils.isStringEmpty(geneRgdIdStr) ) {
                        geneRgdId = Integer.parseInt(geneRgdIdStr);
                    }
                    //gffGeneId = getTokenValue(info, "ID=", ";");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        if( chr==null ) {
                            System.out.println("NULL chr");
                        } else {
                            geneInfo = new GeneInfo();
                            geneInfo.geneSymbol = geneSymbol;
                            geneInfo.ncbiGeneId = ncbiGeneId;
                            geneInfo.geneRgdId = geneRgdId;
                            geneInfo.geneBioType = geneBioType;
                            geneInfo.pseudo = pseudoStr != null && pseudoStr.equals("true");

                            geneInfo.chr = chr;
                            geneInfo.startPos = startPos;
                            geneInfo.stopPos = stopPos;
                            geneInfo.strand = strand;
                            geneMap.put(ncbiGeneId, geneInfo);
                        }
                    } else {
                        throw new Exception("unexpected 1: "+lineNr);
                    }
                }

                case "mRNA", "lnc_RNA", "transcript", "primary_transcript" -> {
                    String ncbiGeneId = getTokenValue(info, "Dbxref=GeneID:", ",");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        throw new Exception("unexpected 2: "+lineNr);
                    }
                    String trAcc = getTokenValue(info, "Name=", ";");
                    String trId = getTokenValue(info, "ID=", ";");
                    // this transcript must be new in trList
                    TrInfo trInfo = null;
                    for( TrInfo ti: geneInfo.trInfos ) {
                        if( ti.id.equals(trId) ) {
                            throw new Exception("unexpected 3: "+lineNr);
                        }
                    }
                    trInfo = new TrInfo();
                    trInfo.id = trId;
                    trInfo.acc = trAcc;
                    trInfo.chr = chr;
                    trInfo.strand = strand;
                    trInfo.startPos = startPos;
                    trInfo.stopPos = stopPos;
                    geneInfo.trInfos.add(trInfo);
                }

                case "exon" -> {
                    String ncbiGeneId = getTokenValue(info, "Dbxref=GeneID:", ",");
                    String trId = getTokenValue(info, "Parent=", ";");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        throw new Exception("unexpected 4: "+lineNr);
                    }
                    TrInfo trInfo = null;
                    for( TrInfo ti: geneInfo.trInfos ) {
                        if( ti.id.equals(trId) ) {
                            ExonInfo exon = new ExonInfo();
                            exon.startPos = startPos;
                            exon.stopPos = stopPos;
                            ti.exons.add(exon);
                            trInfo = ti;
                            break;
                        }
                    }
                    if( trInfo==null ) {
                        if( !ignoredFeatures.contains(trId) ) {
                            // exons without NCBI transcript accessions -- we skip them
                            System.out.println("*** EXON skipped: " + lineNr);
                        }
                    }
                }

                case "CDS" -> {
                    String ncbiGeneId = getTokenValue(info, "Dbxref=GeneID:", ",");
                    String trId = getTokenValue(info, "Parent=", ";");
                    String proteinId = getTokenValue(info, "Name=", ";");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        throw new Exception("unexpected 6: "+lineNr);
                    }
                    TrInfo trInfo = null;
                    for( TrInfo ti: geneInfo.trInfos ) {
                        if( ti.id.equals(trId) ) {
                            trInfo = ti;

                            if( ti.cdsStart==0 ) {
                                ti.cdsStart = startPos;
                            } else if( startPos < ti.cdsStart ) {
                                ti.cdsStart = startPos;
                            }

                            if( ti.cdsStop==0 ) {
                                ti.cdsStop = stopPos;
                            } else if( stopPos > ti.cdsStop ) {
                                ti.cdsStop = stopPos;
                            }

                            if( ti.proteinId==null ) {
                                ti.proteinId = proteinId;
                            }
                            break;
                        }
                    }
                    if( trInfo==null ) {
                        throw new Exception("unexpected 7: "+lineNr);
                    }
                }

                case "miRNA", "antisense_RNA", "snRNA", "rRNA", "telomerase_RNA", "SRP_RNA", "RNase_MRP_RNA" -> {
                    String id = getTokenValue(info, "ID=", ";");
                    ignoredFeatures.add(id);
                }

                case "match", "cDNA_match" -> { // ignore
                }

                default -> {
                    System.out.println("unknown object: " + obj);
                }
            }
        }
        in.close();

        System.out.println("objCount: ");
        for( Map.Entry<String, Integer> entry: objCount.entrySet() ) {
            System.out.println("    "+entry.getKey()+": "+entry.getValue());
        }
        System.out.println("ignored features skipped (miRNA, snRNA, rRNA, antisense_RNA, telomerase_RNA, SRP_RNA): "+ignoredFeatures.size());

        return geneMap;
    }


    void processGene( GeneInfo geneInfo ) throws Exception {

        CounterPool counters = new CounterPool();

        updateGene(geneInfo, counters);

        updateGenePositions(geneInfo, counters);

        updateTranscripts(geneInfo, counters);

        updateTranscriptPositions(geneInfo, counters);

        updateTranscriptsFeatures(geneInfo, counters);

        updateTranscriptVersion(geneInfo, counters);
    }

    void updateGene( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }

    void updateGenePositions( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }

    void updateTranscripts( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }

    void updateTranscriptPositions( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }

    void updateTranscriptsFeatures( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }

    void updateTranscriptVersion( GeneInfo geneInfo, CounterPool counters ) throws Exception {
        throw new Exception("todo");
    }


    String getTokenValue(String info, String startToken, String endToken1, String endToken2) {
        String result = null;
        int p1 = info.indexOf(startToken);
        if( p1>=0 ) {
            int p2 = info.indexOf(endToken1, p1);
            int p3 = info.indexOf(endToken2, p1);
            if( p2>=0 && p3>=0 ) {
                int p4 = Math.min(p2, p3);
                result = info.substring(p1 + startToken.length(), p4);
            }
            else if( p2>=0 && p3<0 ) {
                result = info.substring(p1 + startToken.length(), p2);
            }
            else if( p3>=0 && p2<0 ) {
                result = info.substring(p1 + startToken.length(), p3);
            }
            else {
                result = info.substring(p1 + startToken.length());
            }
        }
        return result;
    }

    String getTokenValue(String info, String startToken, String endToken) {
        String result = null;
        int p1 = info.indexOf(startToken);
        if( p1>=0 ) {
            int p2 = info.indexOf(endToken, p1);
            if( p2>=0 ) {
                result = info.substring(p1 + startToken.length(), p2);
            } else {
                result = info.substring(p1 + startToken.length());
            }
        }
        return result;
    }

    String getTokenValue(String info, String startToken) {
        String result = null;
        int p1 = info.indexOf(startToken);
        if( p1>=0 ) {
            result = info.substring(p1 + startToken.length());
        }
        return result;
    }

    static public class GeneInfo {
        String geneSymbol;
        String ncbiGeneId;
        int geneRgdId;
        String geneBioType;
        boolean pseudo;

        String chr;
        int startPos;
        int stopPos;
        String strand;

        List<TrInfo> trInfos = new ArrayList<>();
    }

    static public class TrInfo {
        String id;  // rna-XM_008004389.1, etc
        String acc; // XM_xxx etc
        String chr;
        String strand; // '+' or '-'
        int startPos;
        int stopPos;

        int cdsStart;
        int cdsStop;
        String proteinId;

        List<ExonInfo> exons = new ArrayList<>();
    }

    static public class ExonInfo {
        int startPos;
        int stopPos;
    }
}
