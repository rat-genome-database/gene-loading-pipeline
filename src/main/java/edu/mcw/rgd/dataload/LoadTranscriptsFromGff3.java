package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class LoadTranscriptsFromGff3 {

    /**
     * usage: LoadTranscriptsFromGff3 [mapKey gff3File]
     * <p>
     * with arguments: loads, or restores, the transcripts of one assembly from an NCBI GFF3 file;
     * f.e. transcripts of mRatBN7.2 from the archived annotation release GCF_015227675.2-RS_2023_06:
     * <pre>LoadTranscriptsFromGff3 372 /data/GCF_015227675.2_mRatBN7.2_genomic.gff.gz</pre>
     * from the pipeline jar, as load_transcripts_from_gff3.sh does:
     * <pre>-jar lib/EntrezGeneLoading.jar -transcripts_from_gff3 372 /data/GCF_...gff.gz -species rat</pre>
     * transcripts already in RGD are matched by accession; transcripts detached in the past are restored
     * under their old rgd id (per STABLE_TRANSCRIPTS); existing feature objects are bound, not duplicated
     * <p>
     * without arguments: loads the strain assemblies listed in run()
     */
    public static void main(String[] args) throws IOException {

        try {
            LoadTranscriptsFromGff3 loader = new LoadTranscriptsFromGff3();
            if( args.length>=2 ) {
                loader.run(Integer.parseInt(args[0]), args[1]);
            } else {
                loader.run();
            }
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    final String SRC_PIPELINE = "NCBI";
    int mapKey;
    CounterPool counters;

    EGDAO dao = EGDAO.getInstance();
    Logger log = LogManager.getLogger("transcripts");

    void run() throws Exception {

        // map_key -> gff file with NCBI pilot annotation, per MAPS.GENBANK_ASSEMBLY_ACC
        Map<Integer, String> gffFiles = new LinkedHashMap<>();
        gffFiles.put(301, "H:/rat_gff3/GCA_023515785.1_UTH_Rnor_SHR_Utx_genomic.gff.gz");
        gffFiles.put(302, "H:/rat_gff3/GCA_021556685.1_UTH_Rnor_SHRSP_BbbUtx_1.0_genomic.gff.gz");
        gffFiles.put(303, "H:/rat_gff3/GCA_023515805.1_UTH_Rnor_WKY_Bbb_1.0_genomic.gff.gz");

        for( Map.Entry<Integer, String> entry: gffFiles.entrySet() ) {
            run(entry.getKey(), entry.getValue());
        }
    }

    void run(int mapKey, String fname) throws Exception {

        this.mapKey = mapKey;
        counters = new CounterPool();

        System.out.println("=== processing mapKey=" + mapKey + ", file " + fname);

        // gene map: ncbi-gene-id -> list of GeneInfo
        Map<String, GeneInfo> geneMap = loadGeneMap(fname);

        List<GeneInfo> randomizedList = new ArrayList<>(geneMap.values());
        Collections.shuffle(randomizedList);

        AtomicInteger i = new AtomicInteger(0);
        AtomicInteger failures = new AtomicInteger(0);
        randomizedList.stream().parallel().forEach( geneInfo -> {

            i.incrementAndGet();
            System.out.println(i+". "+geneInfo.geneSymbol+"   RGD:"+geneInfo.geneRgdId);

            try {
                processGene(geneInfo);
            } catch( Exception e ) {
                // one gene must not abort the whole run: the failure is reported and counted, the run continues
                failures.incrementAndGet();
                counters.increment("GENES: failed with exception");
                log.error("mapKey="+mapKey+" gene "+geneInfo.geneSymbol+" GeneID:"+geneInfo.ncbiGeneId
                        +" RGD:"+geneInfo.geneRgdId+" failed", e);
                e.printStackTrace();
            }

            // dump counters every 1000 genes
            if( i.get()%1000==0 ) {
                System.out.println(counters.dumpAlphabetically());
            }
        });

        System.out.println(counters.dumpAlphabetically());
        if( failures.get()>0 ) {
            System.out.println("WARNING: "+failures.get()+" genes failed with an exception -- see the transcripts log");
        }
    }

    Map<String, GeneInfo> loadGeneMap(String fname) throws Exception {

        BufferedReader in = Utils.openReader(fname);
        String line;
        Map<String, Integer> objCount = new HashMap<>();

        Map<String, GeneInfo> geneMap = new HashMap<>();
        String chr = "", regionChrAcc = "";
        HashSet<String> ignoredFeatures = new HashSet<>(); // we skip these features together with their exon child objects
        HashSet<String> skippedGenes = new HashSet<>(); // NCBI gene ids of genes on unplaced scaffolds (no chromosome) -- skipped with their child features

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
                    // only chromosome-level sequences carry usable coordinates: unplaced scaffolds (genome=genomic)
                    // are labelled with the chromosome they belong to, but their coordinates are scaffold-local
                    String genome = getTokenValue(info, "genome=", ";");
                    chr = genome!=null && genome.equals("genomic") ? null : getTokenValue(info, "Name=", ";");
                    if( chr!=null ) {
                        System.out.println("processing chromosome " + chr);
                    }
                }

                case "pseudogene", "gene" -> {
                    String ncbiGeneId = getGeneId(info);
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
                            // gene on an unplaced scaffold (region without a chromosome name): skip it with its child features
                            skippedGenes.add(ncbiGeneId);
                            counters.increment("GENES: skipped (unplaced scaffold)");
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

                case "mRNA", "lnc_RNA", "transcript", "primary_transcript", "ncRNA", "snoRNA", "snRNA", "rRNA", "tRNA",
                     "miRNA", "antisense_RNA", "telomerase_RNA", "SRP_RNA", "RNase_MRP_RNA", "scRNA", "Y_RNA",
                     "vault_RNA", "guide_RNA" -> {
                    String ncbiGeneId = getGeneId(info);

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        if( skippedGenes.contains(ncbiGeneId) ) {
                            ignoredFeatures.add(getTokenValue(info, "ID=", ";"));
                            break;
                        }
                        throw new Exception("unexpected 2: "+lineNr);
                    }
                    String trAcc = getTokenValue(info, "Name=", ";");
                    String trId = getTokenValue(info, "ID=", ";");

                    // only RefSeq transcripts (NM_, NR_, XM_, XR_) are loaded; records without an accession
                    // (tRNA genes, mature miRNA products) are ignored together with their exons
                    if( trAcc==null || !trAcc.matches("[NX][MR]_[0-9]+(\\.[0-9]+)?") ) {
                        ignoredFeatures.add(trId);
                        counters.increment("TRANSCRIPTS: skipped (no RefSeq accession): "+obj);
                        break;
                    }
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
                    String ncbiGeneId = getGeneId(info);
                    String trId = getTokenValue(info, "Parent=", ";");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        if( skippedGenes.contains(ncbiGeneId) ) {
                            break;
                        }
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
                        if( ignoredFeatures.contains(trId) ) {
                            // exon of an ignored feature (tRNA, mature miRNA, gene segment)
                        } else if( trId!=null && trId.startsWith("gene-") ) {
                            // exon directly under a gene without transcripts (f.e. pseudogene)
                            counters.increment("EXONS: skipped (gene without transcripts)");
                        } else {
                            // exons without NCBI transcript accessions -- we skip them
                            System.out.println("*** EXON skipped: " + lineNr);
                        }
                    }
                }

                case "CDS" -> {
                    String ncbiGeneId = getGeneId(info);
                    String trId = getTokenValue(info, "Parent=", ";");
                    String proteinId = getTokenValue(info, "Name=", ";");

                    GeneInfo geneInfo = geneMap.get(ncbiGeneId);
                    if (geneInfo == null) {
                        if( skippedGenes.contains(ncbiGeneId) ) {
                            break;
                        }
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
                        if( ignoredFeatures.contains(trId) ) {
                            // CDS of an ignored feature (f.e. V_gene_segment)
                            counters.increment("CDS: skipped (ignored parent feature)");
                        } else {
                            throw new Exception("unexpected 7: "+lineNr);
                        }
                    }
                }


                case "V_gene_segment", "C_gene_segment", "D_gene_segment", "J_gene_segment" -> {
                    // immunoglobulin / T-cell receptor gene segments: no transcripts; skipped with their exons and CDS
                    ignoredFeatures.add(getTokenValue(info, "ID=", ";"));
                    counters.increment("GENE SEGMENTS: skipped");
                }

                case "match", "cDNA_match", "D_loop", "origin_of_replication" -> { // ignore
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
        System.out.println("ignored features skipped (tRNA, mature miRNA, gene segments, RNAs without a RefSeq accession): "+ignoredFeatures.size());

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

            updateTranscriptsFeatures(geneInfo, gene);

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

        // no gene matching by EG ID: try to match by the RGD ID given in the gff file, if any
        if( genes.isEmpty() && geneInfo.geneRgdId!=0 ) {
            try {
                genes.add(dao.getGene(geneInfo.geneRgdId));
            } catch( GeneDAO.GeneDAOException e ) {
                // no gene with that rgd id in RGD
            }
        }

        if( genes.isEmpty() ) {
            counters.increment("GENES: no match by EG ID");
            log.info("mapKey="+mapKey+" gene not in RGD, skipped: GeneID:"+geneInfo.ncbiGeneId+" "+geneInfo.geneSymbol
                    +" RGD:"+geneInfo.geneRgdId+" "+geneInfo.geneBioType);
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

        // transcripts are attached to active genes only; a transcript of an inactive gene would be withdrawn by TranscriptQC
        if( gene!=null && !dao.getRgdId(gene.getRgdId()).getObjectStatus().equals("ACTIVE") ) {
            counters.increment("GENES: inactive in RGD, skipped");
            gene = null;
        }

        return gene;
    }

    void updateGenePositions( GeneInfo geneInfo, Gene gene ) throws Exception {

        MapData mdIncoming = new MapData();
        mdIncoming.setMapKey(mapKey);
        mdIncoming.setSrcPipeline(SRC_PIPELINE);
        mdIncoming.setRgdId(gene.getRgdId());
        mdIncoming.setChromosome(geneInfo.chr);
        mdIncoming.setStrand(geneInfo.strand);
        mdIncoming.setStartPos(geneInfo.startPos);
        mdIncoming.setStopPos(geneInfo.stopPos);

        List<MapData> mds = dao.getMapData(gene.getRgdId(), mapKey);
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
                tr.setNonCoding(trInfo.cdsStart==0);

                // accession attached to another gene: not touched, only reported
                List<Transcript> trsByAcc = dao.getTranscriptsByAccId(trAcc);
                if( !trsByAcc.isEmpty() ) {
                    log.warn("mapKey="+mapKey+" "+trAcc+" is attached to gene RGD:"+trsByAcc.get(0).getGeneRgdId()
                            +", not to RGD:"+gene.getRgdId()+" ("+geneInfo.geneSymbol+"): skipped");
                    counters.increment("TRANSCRIPTS: skipped (accession attached to another gene)");
                    trInfo.rgdId = 0;
                    continue;
                }

                // a transcript detached from its gene in the past (row in TRANSCRIPTS deleted, rgd id withdrawn)
                // keeps its rgd id in STABLE_TRANSCRIPTS: that rgd id is reused, so the positions and features
                // it still has on other assemblies are reconnected instead of a new transcript object being created
                int restoredRgdId = restoreWithdrawnTranscript(tr, trAcc);
                if( restoredRgdId!=0 ) {
                    trInfo.rgdId = restoredRgdId;
                    continue;
                }

                dao.createTranscript(tr, SpeciesType.RAT);
                trInfo.rgdId = tr.getRgdId();
                counters.increment("TRANSCRIPTS: inserted");
            }
        }
    }

    /**
     * re-attach a withdrawn transcript: reuse its rgd id from STABLE_TRANSCRIPTS, re-activate the rgd id
     * and insert its row into TRANSCRIPTS
     * @return the reused rgd id, or 0 if the accession has no reusable rgd id
     */
    int restoreWithdrawnTranscript( Transcript tr, String trAcc ) throws Exception {

        for( int rgdId: dao.getTranscriptRgdIdsByAccession(trAcc) ) {
            if( dao.getTranscript(rgdId)!=null ) {
                continue; // rgd id in use by another transcript row
            }
            RgdId id = dao.getRgdId(rgdId);
            if( id==null || id.getObjectKey()!=RgdId.OBJECT_KEY_TRANSCRIPTS || id.getSpeciesTypeKey()!=SpeciesType.RAT ) {
                continue;
            }
            if( !id.getObjectStatus().equals("ACTIVE") ) {
                id.setObjectStatus("ACTIVE");
                id.setLastModifiedDate(new Date());
                dao.updateRgdId(id);
                counters.increment("TRANSCRIPT RGD IDS: re-activated");
            }
            tr.setRgdId(rgdId);
            dao.insertTranscript(tr);
            log.info("mapKey="+mapKey+" restored transcript "+trAcc+" RGD:"+rgdId+" for gene RGD:"+tr.getGeneRgdId());
            counters.increment("TRANSCRIPTS: restored (withdrawn transcript re-attached)");
            return rgdId;
        }
        return 0;
    }

    void updateTranscriptPositions( GeneInfo geneInfo ) throws Exception {

        for( TrInfo trInfo: geneInfo.trInfos ) {
            if( trInfo.rgdId==0 ) {
                continue; // transcript skipped
            }

            MapData mdIncoming = new MapData();
            mdIncoming.setMapKey(mapKey);
            mdIncoming.setSrcPipeline(SRC_PIPELINE);
            mdIncoming.setRgdId(trInfo.rgdId);
            mdIncoming.setChromosome(trInfo.chr);
            mdIncoming.setStrand(trInfo.strand);
            mdIncoming.setStartPos(trInfo.startPos);
            mdIncoming.setStopPos(trInfo.stopPos);

            List<MapData> mds = dao.getMapData(trInfo.rgdId, mapKey);
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

    void updateTranscriptsFeatures( GeneInfo geneInfo, Gene gene ) throws Exception {

        if( geneInfo.geneBioType.equals("miRNA") ) {
            return;
        }

        // feature objects already in RGD for this gene on this assembly, by type and position;
        // an incoming feature matching one of them (f.e. an exon shared with a sibling transcript)
        // is bound to the transcript instead of being created again
        Map<String, Integer> geneFeatureIds = new HashMap<>();
        for( TranscriptFeature f: dao.getFeaturesForGene(gene.getRgdId()) ) {
            if( Integer.valueOf(mapKey).equals(f.getMapKey()) ) {
                geneFeatureIds.putIfAbsent(featureKey(f), f.getRgdId());
            }
        }

        // build incoming features: exons and utrs
        for( TrInfo trInfo: geneInfo.trInfos ) {
            if( trInfo.rgdId==0 ) {
                continue; // transcript skipped
            }

            int trStart = 0;
            int trStop = 0;

            List<TranscriptFeature> features = new ArrayList<>();
            for( ExonInfo exonInfo: trInfo.exons ) {

                MapData md = new MapData();
                md.setMapKey(mapKey);
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
                md.setMapKey(mapKey);
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
                md.setMapKey(mapKey);
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
            List<TranscriptFeature> ftsInRgd = dao.getFeaturesForTr(trInfo.rgdId, mapKey);

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

                if( ftInRgd!=null ) {
                    counters.increment("TR "+ft.getCanonicalName().toUpperCase()+": matched");
                    continue;
                }

                // reuse an existing feature object: one of this gene, or an orphaned one whose transcript
                // was detached in the past (looked up by type and exact position)
                Integer featureRgdId = geneFeatureIds.get(featureKey(ft));
                if( featureRgdId==null ) {
                    List<Integer> ids = dao.getFeatureRgdIdsByPosition(ft);
                    if( !ids.isEmpty() ) {
                        featureRgdId = ids.get(0);
                    }
                }
                if( featureRgdId!=null ) {
                    dao.bindFeatureToTranscript(trInfo.rgdId, featureRgdId);
                    geneFeatureIds.putIfAbsent(featureKey(ft), featureRgdId);
                    counters.increment("TR "+ft.getCanonicalName().toUpperCase()+": bound to existing feature object");
                } else {
                    ft.setTranscriptRgdId(trInfo.rgdId);
                    ft.setSrcPipeline(SRC_PIPELINE);
                    dao.createFeature( ft, SpeciesType.RAT );
                    geneFeatureIds.put(featureKey(ft), ft.getRgdId());
                    counters.increment("TR "+ft.getCanonicalName().toUpperCase()+": inserted");
                }
            }
        }

    }

    String featureKey( TranscriptFeature f ) {
        return f.getFeatureType()+"|"+f.getChromosome()+"|"+f.getStartPos()+"|"+f.getStopPos()+"|"+f.getStrand();
    }

    void updateTranscriptVersion( GeneInfo geneInfo ) throws Exception {

        System.out.println("todo");
    }


    // NCBI gene id from the Dbxref attribute, wherever it is in the Dbxref list (f.e. Dbxref=GenBank:YP_665629.1,GeneID:26193)
    String getGeneId(String info) {
        String dbxref = getTokenValue(info, "Dbxref=", ";");
        return dbxref==null ? null : getTokenValue(dbxref, "GeneID:", ",");
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
