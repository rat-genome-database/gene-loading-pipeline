package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class LoadTranscriptsFromGff3 {

    public static void main(String[] args) throws IOException {

        try {
            new LoadTranscriptsFromGff3().run();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    final int MAP_KEY = 304;
    final String SRC_PIPELINE = "NCBI";
    CounterPool counters = new CounterPool();

    EGDAO dao = EGDAO.getInstance();

    void run() throws Exception {

        String fname = "/Users/mtutaj/Downloads/r/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.gff";
        fname = "/tmp/g/304.gff";

        // gene map: ncbi-gene-id -> list of GeneInfo
        Map<String, GeneInfo> geneMap = loadGeneMap(fname);

        List<GeneInfo> randomizedList = new ArrayList<>(geneMap.values());
        Collections.shuffle(randomizedList);

        AtomicInteger i = new AtomicInteger(0);
        randomizedList.stream().parallel().forEach( geneInfo -> {

            i.incrementAndGet();
            System.out.println(i+". "+geneInfo.geneSymbol+"   RGD:"+geneInfo.geneRgdId);

            try {
                processGene(geneInfo);
            } catch( Exception e ) {
                throw new RuntimeException(e);
            }

            // dump counters every 1000 genes
            if( i.get()%1000==0 ) {
                System.out.println(counters.dumpAlphabetically());
            }
        });

        System.out.println(counters.dumpAlphabetically());
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

                case "antisense_RNA", "miRNA", "snRNA", "rRNA", "telomerase_RNA", "SRP_RNA", "RNase_MRP_RNA" -> {
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

        if( !geneInfo.geneBioType.equals("protein_coding")
         && !geneInfo.geneBioType.equals("pseudogene")
         && !geneInfo.geneBioType.equals("transcribed_pseudogene")
         && !geneInfo.geneBioType.equals("lncRNA")
         && !geneInfo.geneBioType.equals("rRNA")
         && !geneInfo.geneBioType.equals("snRNA")
         && !geneInfo.geneBioType.equals("telomerase_RNA")
         && !geneInfo.geneBioType.equals("misc_RNA")
         && !geneInfo.geneBioType.equals("SRP_RNA")
         && !geneInfo.geneBioType.equals("antisense_RNA")
         && !geneInfo.geneBioType.equals("RNase_MRP_RNA")
         && !geneInfo.geneBioType.equals("miRNA"))
        {
            System.out.println(" not protein coding");
        }

        Gene gene = updateGene(geneInfo);

        if( gene!=null ) {
            updateGenePositions(geneInfo, gene);

            updateTranscripts(geneInfo, gene);

            updateTranscriptPositions(geneInfo);

            updateTranscriptsFeatures(geneInfo);

            //updateTranscriptVersion(geneInfo);
        }
    }

    Gene updateGene( GeneInfo geneInfo ) throws Exception {

        Gene gene = null;
        List<Gene> genes = dao.getGenesByEGID(geneInfo.ncbiGeneId);

        // multi genes: remove inactive genes
        if( genes.size()>1 ) {
            Iterator<Gene> it = genes.iterator();
            while( it.hasNext() ) {
                Gene g = it.next();
                RgdId id = dao.getRgdId(g.getRgdId());
                if( !id.getObjectStatus().equals("ACTIVE") ) {
                    it.remove();
                }
            }
        }
        // multigenes: remove genes with non-matching RGD ID
        if( genes.size()>1 ) {
            genes.removeIf(g -> g.getRgdId() != geneInfo.geneRgdId);
        }

        // no gene matching by EG ID: try to match by RGD ID
        if( genes.isEmpty() ) {
            Gene g = dao.getGene(geneInfo.geneRgdId);
            genes.add(g);
        }

        if( genes.isEmpty() ) {
            counters.increment("GENES: no match by EG ID");
        }
        else if( genes.size()>1 ) {
            counters.increment("GENES: multimatch by EG ID");
        }
        else {
            counters.increment("GENES: single match by EG ID");

            if( genes.get(0).getRgdId() == geneInfo.geneRgdId ) {
                gene = genes.get(0);
                counters.increment("GENES: single match by RGD ID");
            }
            else if( genes.get(0).getSymbol().equals(geneInfo.geneSymbol) ) {
                gene = genes.get(0);
                counters.increment("GENES: single match by symbol");
            }
            else {
                gene = genes.get(0);
                counters.increment("GENES: match by EG ID, but mismatch by RGD ID and symbol");
            }
        }

        return gene;
    }

    void updateGenePositions( GeneInfo geneInfo, Gene gene ) throws Exception {

        MapData mdIncoming = new MapData();
        mdIncoming.setMapKey(MAP_KEY);
        mdIncoming.setSrcPipeline(SRC_PIPELINE);
        mdIncoming.setRgdId(gene.getRgdId());
        mdIncoming.setChromosome(geneInfo.chr);
        mdIncoming.setStrand(geneInfo.strand);
        mdIncoming.setStartPos(geneInfo.startPos);
        mdIncoming.setStopPos(geneInfo.stopPos);

        List<MapData> mds = dao.getMapData(gene.getRgdId(), MAP_KEY);
        for( MapData md: mds ) {
            if( md.equalsByGenomicCoords(mdIncoming) ) {
                counters.increment("GENE POS: matches incoming");
                return;
            }
        }

        // no match: we insert the new pos
        List<MapData> list = new ArrayList<>();
        list.add(mdIncoming);
        if( dao.insertMapData(list) !=0 ) {
            counters.increment("GENE POS: inserted");
        }
    }

    void updateTranscripts( GeneInfo geneInfo, Gene gene ) throws Exception {

        List<Transcript> trsInRgd = dao.getNcbiTranscriptsForGene(gene.getRgdId());

        // find incoming transcript among transcripts in RGD
        for( TrInfo trInfo: geneInfo.trInfos ) {

            boolean trIsInRgd = false;

            int dotPos = trInfo.acc.indexOf(".");
            String trAcc = trInfo.acc.substring(0, dotPos);

            for( Transcript tr: trsInRgd ) {
                if( tr.getAccId().equals(trAcc) ) {
                    counters.increment("TRANSCRIPTS: already in RGD");
                    trInfo.rgdId = tr.getRgdId();
                    trIsInRgd = true;
                    break;
                }
            }
            if( !trIsInRgd ) {
                String proteinAcc = null;
                if( trInfo.proteinId!=null ) {
                    dotPos = trInfo.proteinId.indexOf(".");
                    proteinAcc = trInfo.proteinId.substring(0, dotPos);
                }

                Transcript tr = new Transcript();
                tr.setAccId(trAcc);
                tr.setGeneRgdId(gene.getRgdId());
                tr.setProteinAccId(proteinAcc);
                dao.createTranscript(tr, SpeciesType.RAT);
                trInfo.rgdId = tr.getRgdId();
                counters.increment("TRANSCRIPTS: inserted");
            }
        }
    }

    void updateTranscriptPositions( GeneInfo geneInfo ) throws Exception {

        for( TrInfo trInfo: geneInfo.trInfos ) {

            MapData mdIncoming = new MapData();
            mdIncoming.setMapKey(MAP_KEY);
            mdIncoming.setSrcPipeline(SRC_PIPELINE);
            mdIncoming.setRgdId(trInfo.rgdId);
            mdIncoming.setChromosome(trInfo.chr);
            mdIncoming.setStrand(trInfo.strand);
            mdIncoming.setStartPos(trInfo.startPos);
            mdIncoming.setStopPos(trInfo.stopPos);

            List<MapData> mds = dao.getMapData(trInfo.rgdId, MAP_KEY);
            boolean posIsAlreadyInRgd = false;
            for( MapData md: mds ) {
                if( md.equalsByGenomicCoords(mdIncoming) ) {
                    counters.increment("TR POS: matches incoming");
                    posIsAlreadyInRgd = true;
                    break;
                }
            }

            // no match: we insert the new pos
            if( !posIsAlreadyInRgd ) {
                List<MapData> list = new ArrayList<>();
                list.add(mdIncoming);
                if (dao.insertMapData(list) != 0) {
                    counters.increment("TR POS: inserted");
                }
            }
        }
    }

    void updateTranscriptsFeatures( GeneInfo geneInfo ) throws Exception {

        if( geneInfo.geneBioType.equals("miRNA") ) {
            return;
        }

        // build incoming features: exons and utrs
        for( TrInfo trInfo: geneInfo.trInfos ) {

            int trStart = 0;
            int trStop = 0;

            List<TranscriptFeature> features = new ArrayList<>();
            for( ExonInfo exonInfo: trInfo.exons ) {

                MapData md = new MapData();
                md.setMapKey(MAP_KEY);
                md.setChromosome(trInfo.chr);
                md.setStrand(trInfo.strand);
                md.setStartPos(exonInfo.startPos);
                md.setStopPos(exonInfo.stopPos);

                if( trStart==0 || exonInfo.startPos < trStart ) {
                    trStart = exonInfo.startPos;
                }
                if( trStop==0 || exonInfo.stopPos > trStop ) {
                    trStop = exonInfo.stopPos;
                }

                TranscriptFeature ft = new TranscriptFeature(md);
                ft.setFeatureType(TranscriptFeature.FeatureType.EXON);
                features.add(ft);
            }

            // 1st utr, if available
            if( trInfo.cdsStart!=0 && trStart < trInfo.cdsStart ) {

                MapData md = new MapData();
                md.setMapKey(MAP_KEY);
                md.setChromosome(trInfo.chr);
                md.setStrand(trInfo.strand);
                md.setStartPos(trStart);
                md.setStopPos(trInfo.cdsStart-1);

                TranscriptFeature ft = new TranscriptFeature(md);
                if( trInfo.strand.equals("+") ) {
                    ft.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                } else {
                    ft.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                }
                features.add(ft);
            }

            // 2nd utr, if available
            if( trInfo.cdsStop!=0 && trInfo.cdsStop < trStop ) {

                MapData md = new MapData();
                md.setMapKey(MAP_KEY);
                md.setChromosome(trInfo.chr);
                md.setStrand(trInfo.strand);
                md.setStartPos(trInfo.cdsStop+1);
                md.setStopPos(trStop);

                TranscriptFeature ft = new TranscriptFeature(md);
                if( trInfo.strand.equals("+") ) {
                    ft.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                } else {
                    ft.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                }
                features.add(ft);
            }

            // qc features
            List<TranscriptFeature> ftsInRgd = dao.getFeaturesForTr(trInfo.rgdId, MAP_KEY);

            for( TranscriptFeature ft: features ) {

                // find matching feature in rgd
                TranscriptFeature ftInRgd = null;
                for( TranscriptFeature r: ftsInRgd ) {
                    if (ft.getFeatureType() ==r.getFeatureType()
                    && Utils.intsAreEqual(ft.getStartPos(), r.getStartPos())
                    && Utils.intsAreEqual(ft.getStopPos(), r.getStopPos()) ) {

                        ftInRgd = r;
                        break;
                    }
                }

                if( ftInRgd==null ) {
                    ft.setTranscriptRgdId(trInfo.rgdId);
                    ft.setSrcPipeline(SRC_PIPELINE);
                    dao.createFeature( ft, SpeciesType.RAT );
                    counters.increment("TR "+ft.getCanonicalName().toUpperCase()+": inserted");
                } else {
                    counters.increment("TR "+ft.getCanonicalName().toUpperCase()+": matched");
                }
            }
        }

    }

    void updateTranscriptVersion( GeneInfo geneInfo ) throws Exception {

        System.out.println("todo");
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
        int rgdId;

        List<ExonInfo> exons = new ArrayList<>();
    }

    static public class ExonInfo {
        int startPos;
        int stopPos;
    }
}
