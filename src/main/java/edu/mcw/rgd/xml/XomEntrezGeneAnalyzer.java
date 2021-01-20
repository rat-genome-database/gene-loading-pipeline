package edu.mcw.rgd.xml;

import edu.mcw.rgd.dataload.BulkGene;
import edu.mcw.rgd.dataload.TranscriptInfo;
import edu.mcw.rgd.dataload.TranscriptLocus;
import edu.mcw.rgd.dataload.TranscriptVersionManager;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.*;
import nu.xom.Element;
import org.apache.commons.collections.map.MultiValueMap;
import org.apache.log4j.Logger;
import org.jaxen.JaxenException;
import org.jaxen.XPath;
import org.jaxen.xom.XOMXPath;

import java.io.File;
import java.io.FileWriter;
import java.sql.Timestamp;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since Feb 24, 2010
 * analyze EntrezGene XML file (as downloaded through NCBI eSearch, eFetch utils)
 * using XOM streaming model: parses one record at a time - detect record boundaries ('Entrezgene' elements)
 * then we use XOMXPath objects to extract data useful to us using XPath;
 * to speed-up XPath queries we parse every direct child of 'Entrezgene' element separately
 *
 * to use this class, call method run() passing the File object as parameter;
 * when entire XML record has been parsed
 */
@SuppressWarnings("unchecked")
public class XomEntrezGeneAnalyzer extends XomAnalyzer {

    Logger logAnnotStatus = Logger.getLogger("annot_status");

    // local variables used during parsing xml record
    private BulkGene bulkGene;

    // genomic assemblies to be handled by the pipeline
    java.util.Map<String, String> genomicAssemblies;
    // xpaths associated to the genomic assemblies
    java.util.Map<XPath, Integer> genomicAssemblyXPathMap = new HashMap<>();

    // counts number of occurrences of all genomic assembly names processed by the program
    java.util.Map<String, Integer> genomicAssemblyNameCountMap = new HashMap<>();

    // counts number of occurrences of all assembly names found in the sources
    java.util.Map<String, Integer> anyAssemblyNameCountMap = new HashMap<>();

    // history genomic assemblies (additional positional information provided by RefSeq)
    // maps assembly name to a map key
    java.util.Map<String, Integer> geneLocationHistory;

    // first rec no to be processed
    int firstRecNo;

    // internal record counter
    int recno;

    // name of file to be processed
    String fileName;

    public List<BulkGene> bulkGenes = new ArrayList<>();

    private CounterPool counters;

    // precompiled XPath expressions
    static XPath xpGeneSymbol, xpGeneSymbolAlt, xpGeneName, xpNomenSymbol, xpNomenName;
    static XPath xpGeneOfficialSymbol, xpGeneOfficialName, xpGeneInterimSymbol, xpGeneInterimName, xpGeneSpecies;
    static XPath xpGeneLocHistoryAssembly, xpGeneLocHistoryPos, xpGeneLocHistoryPosMulti, xpGeneLocHistoryChr;
    static XPath xpEntrezgeneID, xpRgdID, xpAliasesS, xpAliasesN, xpAliasesD;
    static XPath xpChromosome, xpChromosome2, xpcM, xpXdbProductAcc;
    static String sAssemblyMapChr, sAssemblyMapScaffold, sAnyAssemblyChr, sAnyAssemblyScaffold;
    static XPath xpAssemblyMapName, xpAssemblyAccession, xpMapStartPos, xpMapStopPos, xpMapStrand, xpAnyAssemblyName;
    static XPath xpAssemblyMapLabel, xpGeneTrackStatus, xpGeneTrackCurrentId, xpXdbMgd, xpXdbHgnc, xpXdbVgnc, xpMTCommentary;
    static XPath xpXdbHomologene, xpXdbKeggReport, xpXdbKeggPathway, xpXdbKeggPathwayName, xpXdbUniProtType;
    static XPath xpXdbUniProt, xpXdbPubMed, xpXdbGenBankN, xpXdbGenBankP, xpXdbEnsembl, xpXdbEnsemblP, xpXdbMiRBase;
    static XPath xpProducts, xpGenomicCoords, xpGenomicCoords2, xpSeqInterval, xpProductType, xpProductLabel, xpAccession, xpAccessionVer;
    static XPath xpProducts2, xpPackedInterval, xpBioSourceGenome;
    static XPath xpGeneRefSeq, xpTrRefSeq, xpTrRefSeqAcc;
    static XPath xpAnnotInfo, xpAnnotLabel, xpAnnotText;

    static {
        try {
            // <Entrezgene_comments>
            xpXdbKeggReport = new XOMXPath(".//Dbtag[Dbtag_db='KEGG']/Dbtag_tag/Object-id/Object-id_str");
            xpXdbKeggPathway = new XOMXPath(".//Dbtag[contains(Dbtag_db,'KEGG pathway')]/Dbtag_tag/Object-id/Object-id_str");
            xpXdbKeggPathwayName = new XOMXPath("../../../../../../../Gene-commentary_text");

            xpXdbUniProt = new XOMXPath(".//Dbtag[starts-with(Dbtag_db,'UniProtKB')]/Dbtag_tag/Object-id/Object-id_str");
            xpXdbUniProtType = new XOMXPath("../../../Dbtag_db");
            xpXdbProductAcc = new XOMXPath("../../../../../../../../../../../Gene-commentary_accession");

            xpXdbPubMed = new XOMXPath(".//Pub/Pub_pmid/PubMedId"); // same query works for _properties section
            xpXdbGenBankN = new XOMXPath(".//Gene-commentary[(contains(Gene-commentary_type/@value,'RNA') or Gene-commentary_type/@value='genomic') and not(starts-with(Gene-commentary_accession,'AC_') or starts-with(Gene-commentary_accession,'NC_') or starts-with(Gene-commentary_accession,'NT_') or starts-with(Gene-commentary_accession,'NW_') or starts-with(Gene-commentary_accession,'NZ_'))]/Gene-commentary_accession");
            xpXdbGenBankP = new XOMXPath(".//Gene-commentary[contains(Gene-commentary_type/@value,'peptide')]/Gene-commentary_accession");
            xpXdbEnsembl = new XOMXPath(".//Dbtag[Dbtag_db='Ensembl']/Dbtag_tag/Object-id/Object-id_str");
            xpGeneRefSeq = new XOMXPath("Gene-commentary[Gene-commentary_heading='RefSeq Status']/Gene-commentary_label");
            xpTrRefSeq = new XOMXPath(".//Gene-commentary/Gene-commentary_comment/Gene-commentary[Gene-commentary_label='RefSeq Status']/Gene-commentary_text");
            xpTrRefSeqAcc = new XOMXPath("../../../Gene-commentary_accession");

            xpAnnotInfo = new XOMXPath("Gene-commentary[Gene-commentary_heading='Annotation Information']/Gene-commentary_properties/Gene-commentary");
            xpAnnotLabel = new XOMXPath("Gene-commentary_label");
            xpAnnotText = new XOMXPath("Gene-commentary_text");

            xpGeneLocHistoryAssembly = new XOMXPath("Gene-commentary[Gene-commentary_heading='Gene Location History']/Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_heading");
            xpGeneLocHistoryChr = new XOMXPath("../Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary");
            xpGeneLocHistoryPos = new XOMXPath("Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval");
            xpGeneLocHistoryPosMulti = new XOMXPath("Gene-commentary_seqs/Seq-loc/Seq-loc_mix/Seq-loc-mix/Seq-loc/Seq-loc_int/Seq-interval");

            // <Entrezgene_gene>
            xpGeneSymbol = new XOMXPath("Gene-ref/Gene-ref_locus");
            xpGeneSymbolAlt = new XOMXPath("Gene-ref/Gene-ref_locus-tag");
            xpGeneName = new XOMXPath("Gene-ref/Gene-ref_desc");
            xpNomenSymbol = new XOMXPath("Gene-ref/Gene-ref_formal-name/Gene-nomenclature/Gene-nomenclature_symbol");
            xpNomenName = new XOMXPath("Gene-ref/Gene-ref_formal-name/Gene-nomenclature/Gene-nomenclature_name");
            xpRgdID = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='RGD']/Dbtag_tag/Object-id/Object-id_id");
            xpAliasesS = new XOMXPath("Gene-ref/Gene-ref_syn/Gene-ref_syn_E");
            xpXdbMgd = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='MGI']/Dbtag_tag/Object-id/Object-id_str");
            xpXdbHgnc = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='HGNC']/Dbtag_tag/Object-id/Object-id_str");
            xpXdbVgnc = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='VGNC']/Dbtag_tag/Object-id/Object-id_str");
            xpXdbEnsemblP = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='Ensembl']/Dbtag_tag/Object-id/Object-id_str");
            xpXdbMiRBase = new XOMXPath("Gene-ref/Gene-ref_db/Dbtag[Dbtag_db='miRBase']/Dbtag_tag/Object-id/Object-id_str");

            // <Entrezgene_homology>
            xpXdbHomologene = new XOMXPath("Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag[contains(Dbtag_db,'HomoloGene')]/Dbtag_tag/Object-id/Object-id_id");

            // <Entrezgene_location>
            xpChromosome = new XOMXPath("Maps[Maps_method/Maps_method_map-type/@value='cyto']/Maps_display-str");
            xpcM = new XOMXPath("Maps[Maps_method/Maps_method_map-type/@value='cM']/Maps_display-str");

            // <Entrezgene_locus>
            sAssemblyMapChr = "Gene-commentary[##ASSEMBLY## and (starts-with(Gene-commentary_accession,'NC_') or starts-with(Gene-commentary_accession,'AC_'))]/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval";
            sAssemblyMapScaffold = "Gene-commentary[##ASSEMBLY## and (contains(Gene-commentary_label,'##SCAFFOLD##'))]/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval";
            xpAssemblyMapName = new XOMXPath("../../../../Gene-commentary_heading");
            xpAssemblyMapLabel = new XOMXPath("../../../../Gene-commentary_label");
            xpAssemblyAccession = new XOMXPath("../Gene-commentary_accession");
            sAnyAssemblyChr = "Gene-commentary[starts-with(Gene-commentary_accession,'NC_') or starts-with(Gene-commentary_accession,'AC_')##SCAFFOLD##]/Gene-commentary_heading";
            sAnyAssemblyScaffold = " or contains(Gene-commentary_label,'##SCAFFOLD##')";
            xpMapStartPos = new XOMXPath("Seq-interval_from");
            xpMapStopPos = new XOMXPath("Seq-interval_to");
            xpMapStrand = new XOMXPath("Seq-interval_strand/Na-strand/@value");
            xpProducts = new XOMXPath("../../../../Gene-commentary_products/Gene-commentary");
            xpGenomicCoords = new XOMXPath("Gene-commentary_genomic-coords/Seq-loc");
            xpGenomicCoords2 = new XOMXPath("Gene-commentary_genomic-coords/Seq-loc/Seq-loc_mix/Seq-loc-mix/Seq-loc");
            xpSeqInterval = new XOMXPath("Seq-loc_int/Seq-interval");
            xpPackedInterval = new XOMXPath("Seq-loc_packed-int/Packed-seqint/Seq-interval");
            xpProductType = new XOMXPath("Gene-commentary_type/@value");
            xpProductLabel = new XOMXPath("Gene-commentary_label");
            xpAccession = new XOMXPath("Gene-commentary_accession");
            xpAccessionVer = new XOMXPath("Gene-commentary_version");
            xpProducts2 = new XOMXPath("Gene-commentary_products/Gene-commentary");

            // mitochondrial gene commentaries do not have headings
            xpMTCommentary = new XOMXPath("Gene-commentary[starts-with(Gene-commentary_accession,'NC_') or starts-with(Gene-commentary_accession,'AC_')]/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval");

            // <Entrezgene_properties>
            xpGeneOfficialSymbol = new XOMXPath("Gene-commentary/Gene-commentary_properties/Gene-commentary[Gene-commentary_label='Official Symbol']/Gene-commentary_text");
            xpGeneOfficialName = new XOMXPath("Gene-commentary/Gene-commentary_properties/Gene-commentary[Gene-commentary_label='Official Full Name']/Gene-commentary_text");
            xpGeneInterimSymbol = new XOMXPath("Gene-commentary/Gene-commentary_properties/Gene-commentary[Gene-commentary_label='Interim Symbol']/Gene-commentary_text");
            xpGeneInterimName = new XOMXPath("Gene-commentary/Gene-commentary_properties/Gene-commentary[Gene-commentary_label='Interim Full Name']/Gene-commentary_text");

            // <Entrezgene_prot>
            xpAliasesN = new XOMXPath("Prot-ref/Prot-ref_name/Prot-ref_name_E");
            xpAliasesD = new XOMXPath("Prot-ref/Prot-ref_desc");

            // <Entrezgene_source>
            xpGeneSpecies = new XOMXPath("BioSource/BioSource_org/Org-ref/Org-ref_taxname");
            xpChromosome2 = new XOMXPath("BioSource/BioSource_subtype/SubSource[SubSource_subtype/@value='chromosome']/SubSource_name");
            xpBioSourceGenome = new XOMXPath("BioSource/BioSource_genome/@value");

            // <Entrezgene_track-info>
            xpEntrezgeneID = new XOMXPath("Gene-track/Gene-track_geneid");
            xpGeneTrackStatus = new XOMXPath("Gene-track/Gene-track_status/@value");
            xpGeneTrackCurrentId = new XOMXPath("Gene-track/Gene-track_current-id/Dbtag[Dbtag_db='GeneID']/Dbtag_tag/Object-id/Object-id_id");
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    //////////////////////////////////////////////////////////
    // SUBRECORD-LEVEL PROCESSING
    //////////////////////////////////////////////////////////

    // analyze <Entrezgene/Entrezgene_comments> element
    void analyzeEntrezgeneComments(Element element) {

        try {

            parseXdb(xpXdbKeggReport, element, XdbId.XDB_KEY_KEGGREPORT);
            parseXdb(xpXdbKeggPathway, element, XdbId.XDB_KEY_KEGGPATHWAY);
            parseXdb(xpXdbUniProt, element, XdbId.XDB_KEY_UNIPROT);
            parseXdb(xpXdbPubMed, element, XdbId.XDB_KEY_PUBMED);
            parseXdb(xpXdbGenBankN, element, XdbId.XDB_KEY_GENEBANKNU);
            parseXdb(xpXdbGenBankP, element, XdbId.XDB_KEY_GENEBANKPROT);
            parseXdbEnsembl(xpXdbEnsembl, element);

            // get gene RefSeq status
            String geneRefSeqStatus = xpGeneRefSeq.stringValueOf(element);
            if( geneRefSeqStatus.length()>0 )
                // set gene refseq status -- make it always uppercase for consistency
                bulkGene.getGene().setRefSeqStatus(geneRefSeqStatus.toUpperCase());

            // get RefSeq status for transcripts
            for( Element el: (List<Element>)xpTrRefSeq.selectNodes(element) ) {

                String refSeqStatus = el.getValue();
                String accId = xpTrRefSeqAcc.stringValueOf(el);

                // set transcript refseq status -- make it always uppercase for consistency
                bulkGene.mapRefSeqStatus.put(accId, refSeqStatus.toUpperCase());
            }

            // annotation information
            for( Element el: (List<Element>)xpAnnotInfo.selectNodes(element) ) {

                String ncbiAnnotStatus = xpAnnotLabel.stringValueOf(el) + ": " + xpAnnotText.stringValueOf(el);
                bulkGene.setNcbiAnnotStatus(ncbiAnnotStatus);
                logAnnotStatus.debug("NCBI_ANNOT_STATUS> EGID:"+this.bulkGene.getEgId()+" "+ncbiAnnotStatus);
            }

            handleGeneLocHistory(element);
        }
        catch(Exception e)
        {
             e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_gene> element
    void analyzeEntrezgeneGene(Element element) {

        try {
            String geneSymbol = xpGeneSymbol.stringValueOf(element);
            if( Utils.isStringEmpty(geneSymbol) ) {
                geneSymbol = xpGeneSymbolAlt.stringValueOf(element);
            }
            bulkGene.getGene().setSymbol(geneSymbol);

            bulkGene.getGene().setName(xpGeneName.stringValueOf(element));
            bulkGene.getGene().setRgdId(xpRgdID.numberValueOf(element).intValue());

            bulkGene.nomenSymbol = xpNomenSymbol.stringValueOf(element);
            bulkGene.nomenName = xpNomenName.stringValueOf(element);

            for( Element el: (List<Element>)xpAliasesS.selectNodes(element) ) {

                Alias alias = new Alias();
                alias.setValue(el.getValue());
                alias.setTypeName("old_gene_symbol");
                alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
                bulkGene.aliases.addAlias(alias);
            }

            parseXdb(xpXdbMgd, element, XdbId.XDB_KEY_MGD);
            parseXdb(xpXdbHgnc, element, XdbId.XDB_KEY_HGNC);
            parseXdb(xpXdbVgnc, element, XdbId.XDB_KEY_VGNC);
            parseXdb(xpXdbMiRBase, element, XdbId.XDB_KEY_MIRBASE);
            parseXdbEnsembl(xpXdbEnsemblP, element);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_homology> element
    void analyzeEntrezgeneHomology(Element element) {

        try {

            parseXdb(xpXdbHomologene, element, XdbId.XDB_KEY_HOMOLOGENE);

        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_locus> element
    void analyzeEntrezgeneLocus(Element element) {

        try {
            // set of encountered assemblies (to avoid duplicate counts when hadnling multiple transcripts for a gene)
            Set<String> encounteredAssemblies = new HashSet<>();

            // mitochondrial genes has simplified parsing
            if( bulkGene.bioSourceGenome.startsWith("mitochondr") ) {
                Element el = (Element) xpMTCommentary.selectSingleNode(element);

                // retrieve the genomic location
                MapData md = new MapData();
                md.setChromosome("MT");
                parseSeqInterval(el, md);
                md.setMapKey(BulkGene.getPrimaryMapKey());
                bulkGene.addMapDataForGene(md);

                bulkGene.setChromosome("MT");
            }
            else
            // run the query for all genomic assemblies and extract map positions
            for( Map.Entry<XPath, Integer> entry: this.genomicAssemblyXPathMap.entrySet() ) {
                XPath xpath = entry.getKey(); 
                List<Element> nodes = xpath.selectNodes(element);
                for( Element el: nodes ) {

                    TranscriptFeature map = createPos(el, TranscriptFeature.FeatureType.TRANSCRIPT, entry.getValue());
                    if( map==null )
                        continue;
                    // add genomic positions
                    bulkGene.addMapDataForGene(map);

                    // check label if it contains chromosome names
                    String chromosomeOverride = checkLabelForChromosome(el, xpAssemblyMapLabel);
                    if( chromosomeOverride!=null )
                        map.setChromosome(chromosomeOverride);

                    // read all gene commentary products
                    List<Element> products = xpProducts.selectNodes(el);
                    for( Element elProduct: products ) {
                        // read transcript accession id
                        String transcriptAccId = xpAccession.stringValueOf(elProduct);
                        if( Utils.isStringEmpty(transcriptAccId) ) {
                            continue;
                        }

                        // transcript accession version
                        String transcriptAccVer = xpAccessionVer.stringValueOf(elProduct);
                        if( !Utils.isStringEmpty(transcriptAccVer) ) {
                            TranscriptVersionManager.getInstance().addVersion(transcriptAccId, transcriptAccVer);
                        }


                        TranscriptInfo transcript = bulkGene.getTranscriptByAccId(transcriptAccId);
                        if( transcript==null )
                            transcript = bulkGene.createTranscript(transcriptAccId);

                        // read product type
                        String productType = xpProductType.stringValueOf(elProduct);
                        if( !productType.equalsIgnoreCase("mRNA")
                         && !productType.equalsIgnoreCase("miscRNA")
                         && !productType.equalsIgnoreCase("pre-RNA")
                         && !productType.equalsIgnoreCase("rRNA")
                         && !productType.equalsIgnoreCase("tRNA")
                         && !productType.equalsIgnoreCase("ncRNA")) {
                            System.out.println("product nonRNA:"+productType+" for "+transcriptAccId);
                            continue;
                        }
                        transcript.setAccId(transcriptAccId);
                        TranscriptLocus locus = transcript.createLocus(map);

                        // read genomic locations
                        // first handle the case when there are multiple CDS regions
                        List<Element> elSeqLocList = (List<Element>) xpGenomicCoords2.selectNodes(elProduct);
                        if( elSeqLocList.size()==0 ) {
                            // try with a single location
                            elSeqLocList = (List<Element>) xpGenomicCoords.selectNodes(elProduct);
                        }

                        // read all genomic locations
                        for( Element elSeqLoc: elSeqLocList ) {
                            Element elSeqInterval = (Element) xpSeqInterval.selectSingleNode(elSeqLoc);
                            if( elSeqInterval!=null ) {
                                TranscriptFeature tfNew = createPos(elSeqInterval, TranscriptFeature.FeatureType.EXON, map.getMapKey());
                                if( chromosomeOverride!=null )
                                    tfNew.setChromosome(chromosomeOverride);
                                locus.addGenomicCoords(tfNew);
                            }
                            for( Element elSeqInterval2: (List<Element>)xpPackedInterval.selectNodes(elSeqLoc) ) {
                                if( elSeqInterval2!=null ) {
                                    TranscriptFeature tfNew = createPos(elSeqInterval2, TranscriptFeature.FeatureType.EXON, map.getMapKey());
                                    if( chromosomeOverride!=null )
                                        tfNew.setChromosome(chromosomeOverride);
                                    locus.addGenomicCoords(tfNew);
                                }
                            }
                        }
                        // set transcript pos be based on the first and last exon
                        // note: transcript could be within a gene region, not covering the entire gene
                        locus.adjustTranscriptPos();

                        // go to mRNA products
                        Element elmRnaProducts = (Element) xpProducts2.selectSingleNode(elProduct);
                        // read mRNA product type
                        String mRnaProductType = xpProductType.stringValueOf(elmRnaProducts);
                        if( elmRnaProducts!=null && mRnaProductType.equalsIgnoreCase("peptide")) {
                            // only mRNA transcripts do have peptide products -- others do not!!!
                            // read peptide accession id
                            String proteinAccId = xpAccession.stringValueOf(elmRnaProducts);
                            if( !proteinAccId.isEmpty() ) {
                                if( transcript.getProteinAccId()!=null && !transcript.getProteinAccId().equals(proteinAccId) )
                                    System.out.println("protein acc id override "+locus.getTranscriptAccId()+" "+locus.getTranscriptCoords()
                                            +" from "+transcript.getProteinAccId()+" to "+proteinAccId);

                                    transcript.setProteinAccId(proteinAccId);
                            }

                            String label = xpProductLabel.stringValueOf(elmRnaProducts);
                            if( !label.isEmpty() ) {
                                counters.increment("PEPTIDE_LABELS");
                                if( transcript.getPeptideLabel()!=null && !transcript.getPeptideLabel().isEmpty() ) {
                                    if( !label.equals(transcript.getPeptideLabel()) )
                                        System.out.println("peptide label override "+locus.getTranscriptAccId()+" "+locus.getTranscriptCoords());
                                }
                                transcript.setPeptideLabel(label);
                            }

                            // read peptide genomic locations
                            // first handle the case when there are multiple CDS regions
                            elSeqLocList = (List<Element>) xpGenomicCoords2.selectNodes(elmRnaProducts);
                            if( elSeqLocList.size()==0 ) {
                                // try with a single location
                                elSeqLocList = (List<Element>) xpGenomicCoords.selectNodes(elmRnaProducts);
                            }

                            // read all peptide genomic locations
                            for( Element elSeqLoc: elSeqLocList ) {
                                Element elSeqInterval = (Element) xpSeqInterval.selectSingleNode(elSeqLoc);
                                if( elSeqInterval!=null ) {
                                    TranscriptFeature tfNew = createPos(elSeqInterval, TranscriptFeature.FeatureType.CDS, map.getMapKey());
                                    if( chromosomeOverride!=null )
                                        tfNew.setChromosome(chromosomeOverride);
                                    locus.addGenomicCoordsForPeptide(tfNew);
                                }
                                for( Element elSeqInterval2: (List<Element>)xpPackedInterval.selectNodes(elSeqLoc) ) {
                                    if( elSeqInterval2!=null ) {
                                        TranscriptFeature tfNew = createPos(elSeqInterval2, TranscriptFeature.FeatureType.CDS, map.getMapKey());
                                        if( chromosomeOverride!=null )
                                            tfNew.setChromosome(chromosomeOverride);
                                        locus.addGenomicCoordsForPeptide(tfNew);
                                    }
                                }
                            }
                        }
                    }

                    // since map positions are extracted, keep count of all assemblies encountered
                    String assemblyName = xpAssemblyMapName.stringValueOf(el);
                    if( !Utils.isStringEmpty(assemblyName) )
                        if( encounteredAssemblies.add(assemblyName) )
                            incrementAssemblyNameCount(assemblyName);
                }
            }

            parseXdb(xpXdbGenBankN, element, XdbId.XDB_KEY_GENEBANKNU);
            parseXdb(xpXdbGenBankP, element, XdbId.XDB_KEY_GENEBANKPROT);

            // find all primary and alternate assembly names found in the sources
            List<Element> nodes = xpAnyAssemblyName.selectNodes(element);
            for( Element el: nodes ) {
                // modify the count of number of appearances of given assembly during the run
                String assemblyName = el.getValue();
                Integer aCount = this.anyAssemblyNameCountMap.get(assemblyName);
                if( aCount==null )
                    aCount = 0;
                this.anyAssemblyNameCountMap.put(assemblyName, 1+aCount);
            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    boolean parseSeqInterval(Element el, MapData md) throws Exception {
        int startPos = 1 + xpMapStartPos.numberValueOf(el).intValue();
        int stopPos = 1 + xpMapStopPos.numberValueOf(el).intValue();
        String strand = xpMapStrand.stringValueOf(el);
        if( strand.equals("plus"))
            strand = "+";
        else if( strand.equals("minus"))
            strand = "-";

        md.setStartPos(startPos);
        md.setStopPos(stopPos);
        md.setStrand(strand);
        return startPos>0 && stopPos>0;
    }

    void handleGeneLocHistory(Element element) throws Exception {

        // multiple assembly patches could be present - they always come from newest to oldest
        // so we keep the most recent position for the assembly
        Set<String> majorAssemblies = new HashSet<>();

        for( Element el: (List<Element>)xpGeneLocHistoryAssembly.selectNodes(element) ) {

            String histAssembly = el.getValue();

            Integer mapKey = this.geneLocationHistory.get(histAssembly);
            if( mapKey==null ) {
                counters.increment("SKIPPED_LOCATIONS_FOR_HISTORIC_ASSEMBLY_"+histAssembly);
            } else if( mapKey!=0 ) {
                // there could be positions for multiple patches of historical assemblies,
                // f.e. "GRCm38", "GRCm38.1", "GRCm38.2"
                // we only load positions from the most recent patch of assembly
                String histMajorAssembly;
                int dotPos = histAssembly.indexOf(".");
                if( dotPos>0 )
                    histMajorAssembly = histAssembly.substring(0, dotPos);
                else
                    histMajorAssembly = histAssembly;

                if( !majorAssemblies.add(histMajorAssembly) )
                    continue; // we already processed this assembly

                // handle primary assembly with multiple patches
                // note: we could already have processed the positions from the latest assembly patch, f.e. GRCm38.p3
                //   and now we could have a different position from lower patch, f.e GRCm38.p2
                //  ensure to not import positions from lower patches!!!
                if( bulkGene.hasGenomicPositionForAssembly(mapKey) ) {
                    continue;
                }

                // sometimes for given assembly, we have both a position on reference assembly
                // and another position on scaffold
                // we load only reference positions
                Set<Integer> scaffoldMapKeys = new HashSet<>(); // no position on ref assembly
                Set<Integer> refMapKeys = new HashSet<>(); // position on ref assembly

                for( Element elChr: (List<Element>) xpGeneLocHistoryChr.selectNodes(el) ) {
                    // accession must start from NC_ or AC_
                    String chrAccId = xpAccession.stringValueOf(elChr);
                    if( !chrAccId.startsWith("NC_") && !chrAccId.startsWith("AC_") ) {
                        scaffoldMapKeys.add(mapKey);
                        continue;
                    }
                    else
                        refMapKeys.add(mapKey);

                    String chr = checkLabelForChromosome(elChr, xpAnnotLabel);

                    // handle single locus
                    handleHistoryPos(xpGeneLocHistoryPos, elChr, chr, mapKey);

                    // handle multiple loci
                    handleHistoryPos(xpGeneLocHistoryPosMulti, elChr, chr, mapKey);
                }

                for( Integer scaffoldMapKey: scaffoldMapKeys ) {
                    if( refMapKeys.contains(scaffoldMapKey) ) {
                        // there was both a scaffold and ref key position for given map
                        // since there was a ref map pos, everything was good
                        continue;
                    }

                    // there was no position on ref map key, only scaffold pos
                    // if we have a position in RGD for this map key, this will mean inconsistency
                    if( bulkGene.getDao().getMapDataCount(bulkGene.getEgId(), mapKey)>0 ) {
                        FileWriter writer = new FileWriter("data/histLocBad.txt", true);
                        writer.write("EGID="+this.bulkGene.getEgId()+",MAP_KEY="+mapKey+"\n");
                        writer.close();
                    }
                }
            }
        }
    }

    void handleHistoryPos(XPath xpGeneLocPos, Element elChr, String chr, Integer mapKey) throws Exception {
        for( Element elLocus: (List<Element>) xpGeneLocPos.selectNodes(elChr) ) {
            MapData md = new MapData();
            md.setChromosome(chr);
            parseSeqInterval(elLocus, md);
            md.setMapKey(mapKey);
            bulkGene.addMapDataForGene(md);
        }
    }

    // analyze <Entrezgene/Entrezgene_location> element
    void analyzeEntrezgeneLocation(Element element) {

        try {
            // extract the cyto string, which looks like '10q24' or 'Xcen-q11'
            // the first part, is a 'chromosome' (like '10' or 'X')
            // the rest is a 'fishband' (like 'q24' or 'cen-q11')
            String raw = xpChromosome.stringValueOf(element);
            String chromosome = "";
            int i;
            for( i=0; i<raw.length(); i++ ) {
                char c = raw.charAt(i);
                // chromosome character must be a digit or one of uppercase letters: 'X','Y','M','T'
                if( Character.isDigit(c) || c=='X' || c =='Y' || c=='M' || c=='T' )
                    chromosome += c;
                else
                    break;
            }
            bulkGene.setChromosome(chromosome);
            bulkGene.setFishband(raw.substring(i).trim());

            // try cM map
            raw = xpcM.stringValueOf(element);
            if( raw.length()>0 ) {
                // we should get something like: "4 37.8 cM"
                int at1 = raw.indexOf(' ');
                int at2 = raw.lastIndexOf(' ');
                if( at1>0 && at2>0 && at2>at1 ) {
                    bulkGene.setChromosome(raw.substring(0, at1));
                    bulkGene.setcM(Double.parseDouble(raw.substring(at1+1, at2)));
                }
            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_properties> element
    void analyzeEntrezgeneProperties(Element element) {

        try {
            bulkGene.geneOfficialSymbol = xpGeneOfficialSymbol.stringValueOf(element);
            bulkGene.geneOfficialName = xpGeneOfficialName.stringValueOf(element);

            bulkGene.geneInterimSymbol = xpGeneInterimSymbol.stringValueOf(element);
            bulkGene.geneInterimName = xpGeneInterimName.stringValueOf(element);

            parseXdb(xpXdbPubMed, element, XdbId.XDB_KEY_PUBMED);
            parseXdb(xpXdbGenBankN, element, XdbId.XDB_KEY_GENEBANKNU);
            parseXdb(xpXdbGenBankP, element, XdbId.XDB_KEY_GENEBANKPROT);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_prot> element
    void analyzeEntrezgeneProt(Element element) {

        try {
            // 'old_gene_name' aliases
            for( Element el: (List<Element>)xpAliasesN.selectNodes(element) ) {

                Alias alias = new Alias();
                alias.setValue(el.getValue());
                alias.setTypeName("old_gene_name");
                alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
                bulkGene.aliases.addAlias(alias);
            }

            // handle 'Preferred protein name', in case it differs from one of the above
            String aliasD = xpAliasesD.stringValueOf(element);
            if( !Utils.isStringEmpty(aliasD) ) {

                // try to add it to aliases (duplicates are not added)
                Alias alias = new Alias();
                alias.setValue(aliasD);
                alias.setTypeName("old_gene_name");
                alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
                bulkGene.aliases.addAlias(alias);
            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_source> element
    void analyzeEntrezgeneSource(Element element) {

        try {
            bulkGene.getGene().setSpeciesTypeKey(SpeciesType.parse(xpGeneSpecies.stringValueOf(element)));
            bulkGene.chromosome2 = xpChromosome2.stringValueOf(element);
            bulkGene.bioSourceGenome = xpBioSourceGenome.stringValueOf(element);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_summary> element
    void analyzeEntrezgeneSummary(Element element) {

        bulkGene.getGene().setDescription(element.getValue().trim());
    }

    // analyze <Entrezgene/Entrezgene_track-info> element
    void analyzeEntrezgeneTrackInfo(Element element) {

        try {
            bulkGene.setEgId(xpEntrezgeneID.stringValueOf(element));

            String geneTrackStatus = xpGeneTrackStatus.stringValueOf(element).toUpperCase();
            if( geneTrackStatus.isEmpty() )
                geneTrackStatus = "UNKNOWN";
            bulkGene.setGeneTrackStatus(geneTrackStatus);

            String geneTrackCurrentId = xpGeneTrackCurrentId.stringValueOf(element);
            bulkGene.setGeneTrackCurrentId(geneTrackCurrentId);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    // analyze <Entrezgene/Entrezgene_type> element
    void analyzeEntrezgeneType(Element element) {

        bulkGene.setEgType(element.getAttributeValue("value"));
    }

    // helper function used to extract a variable list of XDB keys of particular type
    void parseXdb(XPath xpath, Element element, int xdbKey) throws Exception {
        List<Element> nodes = xpath.selectNodes(element);
        for( Element el : nodes ) {
            XdbId xdb = new XdbId();
            xdb.setAccId(el.getValue());
            xdb.setXdbKey(xdbKey);

            // extract KEGG pathway name if applicable
            if( xdbKey==XdbId.XDB_KEY_KEGGPATHWAY ) {
                String pathwayName = xpXdbKeggPathwayName.stringValueOf(el);
                if( pathwayName.startsWith("KEGG pathway: ") )
                    xdb.setLinkText(pathwayName.substring(14));
            }
            // extract Uniprot db name, Trembl or SwissProt, and put it into notes
            else if( xdbKey==XdbId.XDB_KEY_UNIPROT ) {
                // ensure accession id does not have a version suffix
                // f.e. convert 'Q9CQV8.3' into 'Q9CQV8'
                int dotPos = xdb.getAccId().lastIndexOf('.');
                if( dotPos>=0 ) {
                    xdb.setAccId(xdb.getAccId().substring(0, dotPos));
                }

                String uniProtDbType = xpXdbUniProtType.stringValueOf(el);
                xdb.setNotes(uniProtDbType);
                String productAccId = xpXdbProductAcc.stringValueOf(el);
                if( !Utils.isStringEmpty(productAccId) ) {
                    bulkGene.transcriptXdbIds.addUniProtId(xdb.getAccId(), uniProtDbType, productAccId);
                }
            }
            bulkGene.addXdbId(xdb);
        }
    }

    // helper function used to extract Ensembl IDs
    void parseXdbEnsembl(XPath xpath, Element element) throws Exception {
        List<Element> nodes = xpath.selectNodes(element);
        for( Element el : nodes ) {

            String accId = el.getValue();
            // ACC_ID must be not null
            if( accId!=null && accId.length()>7 && accId.startsWith("ENS") ) {
                // determine particular kind of ENSEMBL ID, for rat, mouse and human
                int xdbKey = 0;
                if( accId.startsWith("ENSG0") || accId.charAt(6)=='G') {
                    xdbKey = XdbId.XDB_KEY_ENSEMBL_GENES;
                } else
                if( accId.startsWith("ENSP0") || accId.charAt(6)=='P') {
                    xdbKey = XdbId.XDB_KEY_ENSEMBL_PROTEIN;
                } else
                if( accId.startsWith("ENST0") || accId.charAt(6)=='T') {
                    xdbKey = XdbId.XDB_KEY_ENSEMBL_TRANSCRIPT;
                }

                // add it
                if( xdbKey>0 ) {

                    // remove version from acc id; f.e. 'ENSRNOT00000025450.7' -> 'ENSRNOT00000025450'
                    int dotPos = accId.lastIndexOf('.');
                    if( dotPos>0 ) {
                        accId = accId.substring(0, dotPos);
                    }

                    XdbId xdb = new XdbId();
                    xdb.setAccId(accId);
                    xdb.setXdbKey(xdbKey);
                    bulkGene.addXdbId(xdb);
                }
            }
        }
    }

    TranscriptFeature createPos(Element el, TranscriptFeature.FeatureType feature, int mapKey) throws Exception {

        String notes = "From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString();
        TranscriptFeature fd = new TranscriptFeature();
        if( !parseSeqInterval(el, fd) )
            return null;
        fd.setMapKey(mapKey);
        fd.setNotes(notes);
        fd.setFeatureType(feature);
        return fd;
    }
    //////////////////////////////////////////////////////////
    // RECORD-LEVEL PROCESSING
    //////////////////////////////////////////////////////////

    // overridden to keep track of current element depth
    public void initRecord(String name) {

        // just start new bulk gene
        bulkGene = new BulkGene();
        bulkGene.setRecNo(++recno); // set unique record number
    }

    // entire <EntrezGene> record had been parsed; prepare data for core analysis
    public Element parseRecord(Element element) {

        if( !skipCurrentRecord() ) {
            Gene gene = bulkGene.getGene();

            // alternate chromosome source (to be used for orthologs)
            if( bulkGene.getChromosome()==null && !Utils.isStringEmpty(bulkGene.chromosome2) )
                bulkGene.setChromosome(bulkGene.chromosome2);

            // create cytogenetic mapping for a chromosome and mandatory fishband
            if( !Utils.isStringEmpty(bulkGene.getChromosome()) && !Utils.isStringEmpty(bulkGene.getFishband())
                    && genomicAssemblies!=null && genomicAssemblies.get("Cytomap")!=null ) {

                MapData map = new MapData();
                map.setChromosome(bulkGene.getChromosome());
                map.setFishBand(bulkGene.getFishband());
                map.setNotes("Created on " + (new Timestamp(System.currentTimeMillis())).toString());
                map.setMapKey(Integer.parseInt(genomicAssemblies.get("Cytomap")));
                bulkGene.genePositions.addMapData(map);
            }

            // create cM mapping for a chromosome and optional cM value
            if( !Utils.isStringEmpty(bulkGene.getChromosome()) && bulkGene.getcM()!=null ) {
                MapData map = new MapData();
                map.setChromosome(bulkGene.getChromosome());
                map.setAbsPosition(bulkGene.getcM());
                map.setNotes("Created on " + (new Timestamp(System.currentTimeMillis())).toString());
                String assemblyMapKey = genomicAssemblies.get("cM");
                if( !Utils.isStringEmpty(assemblyMapKey) ) {
                    map.setMapKey(Integer.parseInt(assemblyMapKey));
                    bulkGene.genePositions.addMapData(map);
                }
            }

            // if chromosome information is missing (possible with real data)
            // set the chromosome to "Un"
            if( Utils.isStringEmpty(bulkGene.getChromosome()))
                bulkGene.setChromosome("Un");

            // create alias equal to contig accession id for scaffold genes
            if( bulkGene.getChromosome().length()> 6 ) {

                // the alias must be different than 'old_gene_symbol' or 'old_gene_name' so it won't be shown on gene report pages
                Alias alias = new Alias();
                alias.setTypeName("alternate_id");
                alias.setValue(bulkGene.getChromosome());
                bulkGene.aliases.addAlias(alias);
            }

            // update mappings for chromosome
            bulkGene.genePositions.updateChromosome(bulkGene.getChromosome());

            // update chromosome and refseq_status for transcripts
            for( TranscriptInfo ti: bulkGene.transcripts ) {
                ti.updateChromosome(bulkGene.getChromosome());
                String refSeqStatus = bulkGene.mapRefSeqStatus.get(ti.getAccId());
                if( refSeqStatus!=null )
                    ti.setRefSeqStatus(refSeqStatus);
            }

            // set gene type
            if( bulkGene.getEgType() != null && bulkGene.getEgType().length()>0 )
                gene.setType(bulkGene.getEgType());

            determineGeneSymbol(gene);
            determineGeneName(gene);

            // set notes
            gene.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());

            // add entrez gene id to the list of XDB IDS updated/inserted by pipeline
            XdbId eg = new XdbId();
            eg.setXdbKey(XdbId.XDB_KEY_ENTREZGENE);
            eg.setAccId(bulkGene.getEgId());
            eg.setLinkText(gene.getSymbol());
            bulkGene.addXdbId(eg);

            // xdb ids -- compute linkText and srcPipeline
            for( XdbId xdb : bulkGene.getXdbIds() ) {
                xdb.setSrcPipeline(XdbManager.EG_PIPELINE);
            }

            // mappings done, we can create XML representation of the gene
            bulkGene.setXmlString(bulkGene.toXmlString());

            // add the current gene name and symbol as aliases
            // alias QC will remove any redundancies
            Alias alias = new Alias();
            alias.setTypeName("old_gene_name");
            alias.setValue(gene.getName());
            bulkGene.aliases.addAlias(alias);

            alias = new Alias();
            alias.setTypeName("old_gene_symbol");
            alias.setValue(gene.getSymbol());
            bulkGene.aliases.addAlias(alias);

            //System.out.println("XML recno "+bulkGene.getRecNo());
        }
        else {
            // no XML representation
            bulkGene.setXmlString("");

            System.out.println("XML record skipped: "+bulkGene.getRecNo());
        }

        // add BulkGene to queue for further processing
        bulkGenes.add(bulkGene);

        return null; // discard the element from resulting document
    }

    void determineGeneSymbol(Gene gene) {

        // nomenclature symbols always are preferred
        if( !Utils.isStringEmpty(bulkGene.nomenSymbol) ) {
            gene.setSymbol(bulkGene.nomenSymbol);
            return;
        }

        // if nomenclature symbol is not available, use the preferred symbol
        if( !Utils.isStringEmpty(bulkGene.geneOfficialSymbol) ) {
            gene.setSymbol(bulkGene.geneOfficialSymbol);
            return;
        }

        // if preferred symbol is not available, use the interim symbol
        if( !Utils.isStringEmpty(bulkGene.geneInterimSymbol) ) {
            gene.setSymbol(bulkGene.geneInterimSymbol);
            return;
        }

        // if the gene symbol is still empty, use the first alias of type 'old_gene_symbol' as gene symbol
        if( Utils.isStringEmpty(gene.getSymbol()) && bulkGene.aliases.getIncoming().size() > 0) {

            // the alias must be of 'old_gene_symbol' type to preserve compatibility with old version
            Alias alias = bulkGene.aliases.getIncoming().get(0);
            if( alias.getTypeName().equals("old_gene_symbol") ) {
                gene.setSymbol(alias.getValue());
                bulkGene.aliases.getIncoming().remove(0);
            }
        }

        if( Utils.isStringEmpty(gene.getSymbol())) {
            System.out.println("empty gene symbol");
        }
    }

    void determineGeneName(Gene gene) {

        // gene nomenclature name is always preferred
        if( !Utils.isStringEmpty(bulkGene.nomenName) ) {
            gene.setName(bulkGene.nomenName);
            return;
        }

        // if nomenclature name is not available, use the preferred name
        if( !Utils.isStringEmpty(bulkGene.geneOfficialName) ) {
            gene.setName(bulkGene.geneOfficialName);
            return;
        }

        // if preferred name is not available, use the interim name
        if( !Utils.isStringEmpty(bulkGene.geneInterimName) ) {
            gene.setName(bulkGene.geneInterimName);
            return;
        }

        // if the gene name is still empty, use the first alias of type 'old_gene_name' as gene name
        if( Utils.isStringEmpty(gene.getName()) && bulkGene.aliases.getIncoming().size() > 0) {

            // the alias must be of 'old_gene_name' type to preserve compatibility with old version
            Alias alias = bulkGene.aliases.getIncoming().get(0);
            if( alias.getTypeName().equals("old_gene_name") ) {
                gene.setName(alias.getValue());
                bulkGene.aliases.getIncoming().remove(0);
            }
        }
    }

    public Element parseSubrecord(Element element) {

        if( !skipCurrentRecord() ) {
            switch (element.getLocalName()) {
                case "Entrezgene_comments":
                    analyzeEntrezgeneComments(element);
                    break;
                case "Entrezgene_gene":
                    analyzeEntrezgeneGene(element);
                    break;
                case "Entrezgene_homology":
                    analyzeEntrezgeneHomology(element);
                    break;
                case "Entrezgene_location":
                    analyzeEntrezgeneLocation(element);
                    break;
                case "Entrezgene_locus":
                    analyzeEntrezgeneLocus(element);
                    break;
                case "Entrezgene_properties":
                    analyzeEntrezgeneProperties(element);
                    break;
                case "Entrezgene_prot":
                    analyzeEntrezgeneProt(element);
                    break;
                case "Entrezgene_source":
                    analyzeEntrezgeneSource(element);
                    break;
                case "Entrezgene_summary":
                    analyzeEntrezgeneSummary(element);
                    break;
                case "Entrezgene_track-info":
                    analyzeEntrezgeneTrackInfo(element);
                    break;
                case "Entrezgene_type":
                    analyzeEntrezgeneType(element);
                    break;
            }
        }
        return null;
    }

    // return true if the current record is to be skipped
    boolean skipCurrentRecord() {
        return bulkGene.getRecNo() < this.getFirstRecNo();
    }

    //////////////////////////////////////////////////////////
    // GETTERS - SETTERS
    //////////////////////////////////////////////////////////

    public java.util.Map<String, String> getGenomicAssemblies() {
        return genomicAssemblies;
    }

    public void setGenomicAssemblies(java.util.Map<String, String> genomicAssemblies, String scaffoldAssembly) throws JaxenException {
        this.genomicAssemblies = genomicAssemblies;

        // reverse genomic assemblies map, i.e. build map of map-keys mapped to one or more assembly names
        MultiValueMap rmap = new MultiValueMap();
        if( genomicAssemblies!=null ) {
            for (java.util.Map.Entry<String, String> entry : genomicAssemblies.entrySet()) {
                rmap.put(entry.getValue(), entry.getKey());
            }
        }

        String assembly, anyAssembly;
        if( scaffoldAssembly!=null ) {
            assembly = sAssemblyMapScaffold.replace("##SCAFFOLD##", scaffoldAssembly);
            anyAssembly = sAnyAssemblyScaffold.replace("##SCAFFOLD##", scaffoldAssembly);
        } else {
            assembly = sAssemblyMapChr;
            anyAssembly = "";
        }
        anyAssembly = sAnyAssemblyChr.replace("##SCAFFOLD##", anyAssembly);
        xpAnyAssemblyName = new XOMXPath(anyAssembly);

        // traverse all maps and build XPATH for matching multiple assembly map keywords
        for( Object mapKey: rmap.keySet() ) {
            List<String> assemblyNames = (List<String>) rmap.get(mapKey);
            String xpathExp = "";
            for( String assemblyName: assemblyNames ) {
                if( !xpathExp.isEmpty() )
                    xpathExp += " or ";
                xpathExp += "(contains(Gene-commentary_heading,'"+assemblyName+"')";
                xpathExp += "  and string-length(Gene-commentary_heading)="+assemblyName.length()+") ";
            }
            xpathExp = "("+xpathExp+")";

            genomicAssemblyXPathMap.put(new XOMXPath(assembly.replace("##ASSEMBLY##", xpathExp)), Integer.parseInt(mapKey.toString()));
        }
    }

    void incrementAssemblyNameCount(String assemblyName) {
        Integer oldCount = genomicAssemblyNameCountMap.get(assemblyName);
        if( oldCount == null )
            oldCount = 0;
        genomicAssemblyNameCountMap.put(assemblyName, oldCount+1);
    }

    public Map<String, Integer> getGenomicAssemblyNameCountMap() {
        return genomicAssemblyNameCountMap;
    }

    public void setGenomicAssemblyNameCountMap(Map<String, Integer> genomicAssemblyNameCountMap) {
        this.genomicAssemblyNameCountMap = genomicAssemblyNameCountMap;
    }

    public int getFirstRecNo() {
        return firstRecNo;
    }

    public void setFirstRecNo(int firstRecNo) {
        this.firstRecNo = firstRecNo;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public Map<String, Integer> getAnyAssemblyNameCountMap() {
        return anyAssemblyNameCountMap;
    }

    public void setAnyAssemblyNameCountMap(Map<String, Integer> anyAssemblyNameCountMap) {
        this.anyAssemblyNameCountMap = anyAssemblyNameCountMap;
    }

    String checkLabelForChromosome(Element el, XPath xpath) throws JaxenException {

        String label = xpath.stringValueOf(el);
        String[] words = label.split("\\s+");
        // go over words until we find one named chromosome
        for( int i=0; i<words.length; i++ ) {
            if( words[i].equalsIgnoreCase("chromosome") ) {
                // found 'chromosome' word -- any word after it should be a chromosome name
                if( i<words.length-1 ) {
                    String chr = words[i+1].toUpperCase();
                    // accept any 1..2-char long chromosome
                    if( chr.length()>0 && chr.length()<=2 ) {
                        return chr; // valid chromosome
                    }
                }
            }
        }

        // special handling for scaffold-only assemblies
        String lcScaffoldPattern = " unplaced scaffold reference";
        int scaffoldPatternPos = label.toLowerCase().indexOf(lcScaffoldPattern);
        if( scaffoldPatternPos > 0 ) {
            String chr = label.substring(0, scaffoldPatternPos);
            // for scaffold acc ids, do *NOT* load the version nr,
            // f.e. 'NW_004936963' is fine, 'NW_004936963.1' is not
            if( chr.contains(".") ) {
                chr = chr.substring(0, chr.lastIndexOf('.'));
            }
            this.bulkGene.setChromosome(chr);
            return chr;
        }

        return null;// not found chromosome in label
    }

    public void setGeneLocationHistory(Map<String, Integer> geneLocationHistory) {
        this.geneLocationHistory = geneLocationHistory;
    }

    public CounterPool getCounters() {
        return counters;
    }

    public void setCounters(CounterPool counters) {
        this.counters = counters;
    }


    public List<BulkGene> process(CounterPool counters) throws Exception {

        setCounters(counters);

        // use non-validating xml parser
        setValidate(false);

        parse(new File(getFileName()));

        return bulkGenes;
    }
}
