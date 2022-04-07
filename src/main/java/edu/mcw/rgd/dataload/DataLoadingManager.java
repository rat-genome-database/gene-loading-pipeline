package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.log.*;
import edu.mcw.rgd.process.*;
import edu.mcw.rgd.process.mapping.MapManager;
import edu.mcw.rgd.xml.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map;

/**
 * @author dli, revamped by Marek Tutaj
 *
 * Loading Entrezgene into RGD database
 * The DataLoadingManager contains main method to run the whole pipeline in different user options.
 * 
 */
public class DataLoadingManager {

    private EntrezGeneExtractor entrezGeneExtractor;
    private QualityCheckBulkGene qualityCheck;
    private DecisionMaker decisionMaker;
    private RGDSpringLogger rgdLogger;
    private NcbiEutils eUtils;
    // collection of eg ids processed by the pipeline
    private Map<Integer,Integer> egIds = new HashMap<>();
    private int speciesTypeKey = SpeciesType.ALL; // species type unknown

    private GeneStatusIssueTracker geneStatusIssueTracker = new GeneStatusIssueTracker();

    private Map<String,Map<String,String>> genomicAssemblies;
    private Map<String,String> geneLocationHistory;
    int firstRecNo; // number of record from which the data should be processed
    int qualityCheckingThreadCount; // number of quality checking threads to be run in parallel
    int bgTotal;
    int genesWithWrongType=0;

    PipelineLogger dbLogger = PipelineLogger.getInstance();
    PipelineLogFlagManager dbFlagManager = new PipelineLogFlagManager(dbLogger);

    static long startMilisec=System.currentTimeMillis();
    static long runSec=0;
    
    protected final Logger logger = LogManager.getLogger("process");
    private String version;
    private Map<Integer, String> scaffoldAssemblies;

    public static DataLoadingManager getInstance() {
        return _instance;
    }
    private static DataLoadingManager _instance;

	public static void main(String args[]) {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        DataLoadingManager manager =(DataLoadingManager) (bf.getBean("dataLoadingManager"));
        // write out the current version and connection information
        System.out.println(manager.getVersion());
        System.out.println(manager.rgdLogger.getConnectionInfo());

        if (args.length>=1) {
            try {

            // run the entrezgene file, pass the entrezgene file name
            if (args[0].contains("process_only") && args.length>=4 ){
                manager.initDbLog(manager.getSpecies(args, 2), "process_only", args[1]);

                manager.runEntrezGeneFile(args[1]);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
                manager.geneStatusIssueTracker.writeIssuesToFile();
            }
            // download the entrezgene data, pass the date
            else if (args[0].contains("download_only") && args.length>=5) {
                manager.initDbLog(manager.getSpecies(args, 3), "download_only", args[1]+"-"+args[2]);

                manager.download(args[1], args[2]);
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            // download the data and then load the data
            else if (args[0].contains("download+process") && args.length>=5) {
                manager.initDbLog(manager.getSpecies(args, 3), "download+process", args[1]+"-"+args[2]);

                manager.run(args[1], args[2]);
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();                
            }
            // download the mitochondrial genes and then load the data
            else if (args[0].contains("mitochondrial") && args.length>=3) {
                manager.entrezGeneExtractor.setMitochondrialGenesLoad(true);
                manager.initDbLog(manager.getSpecies(args, 1), "mitochondrial", args[2]);

                manager.run("auto", "auto");
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
            }
            // download and process microRNA genes
            else if (args[0].contains("microRNA") && args.length>=3) {
                manager.initDbLog(manager.getSpecies(args, 1), "microRNA", args[2]);

                String fileName = manager.entrezGeneExtractor.downloadAndProcessGeneList("microRNA", 
                        EGDAO.getInstance().getEdIdsForMicroRnaGenes(manager.speciesTypeKey));
                manager.parseEntrezGeneFile(fileName, null);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
            }
            // download and process microRNA genes
            else if (args[0].contains("tRNA") && args.length>=3) {
                manager.initDbLog(manager.getSpecies(args, 1), "tRNA", args[2]);

                String fileName = manager.entrezGeneExtractor.downloadAndProcessGeneList("tRNA",
                        EGDAO.getInstance().getEdIdsForTRnaGenes(manager.speciesTypeKey));
                manager.parseEntrezGeneFile(fileName, null);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
            }
            else if (args[0].contains("merged") && args.length>=3) {
                manager.initDbLog(manager.getSpecies(args, 1), "merged", args[2]);

                String fileName = manager.entrezGeneExtractor.downloadAndProcessGeneList("merged",
                        EGDAO.getInstance().getEdIdsForMergedGenes());
                manager.parseEntrezGeneFile(fileName, null);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
            }
            // download the mitochondrial genes and then load the data
            else if (args[0].contains("no_strand") && args.length>=4) {

                manager.initDbLog(manager.getSpecies(args, 2), "no_strand", args[3]);

                int mapKey = Integer.parseInt(args[1]);

                manager.entrezGeneExtractor.eSearchOverride = new NcbiEutils().createESearchResult();
                manager.entrezGeneExtractor.eSearchOverride.ids =
                    EGDAO.getInstance().getEgIdsForGenesWithoutStrand(manager.speciesTypeKey, mapKey);
                manager.entrezGeneExtractor.eSearchOverride.recordCount = Integer.toString(manager.entrezGeneExtractor.eSearchOverride.ids.size());
                manager.entrezGeneExtractor.eSearchOverride.queryTranslation = "_no_strand_species_"+manager.speciesTypeKey+"_map_key_"+mapKey+"_";

                manager.run("auto", "auto");
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
            }
            // transform NCBI XML file into EG XML file (BulkGene XML file)
            else if (args[0].contains("transform") && args.length>=5) {
                manager.initDbLog(manager.getSpecies(args, 3), "transform", args[1]+"-"+args[2]);
                manager.getSpecies(args, 3);

                manager.transform(args[1], args[2]);
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            // compare list of removed RefSeq acc ids against RGD
            else if (args[0].contains("refseq_removed")) {

                RefSeqRemoved validator = (RefSeqRemoved) (bf.getBean("refSeqRemoved"));
                validator.run();
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            //
            else if (args[0].contains("refseq_qc_protein")) {

                RefSeqQcProtein qc = new RefSeqQcProtein();
                int mapKey = 360;
                if( args.length>=2 && args[1].startsWith("-mapKey=") ) {
                    mapKey = Integer.parseInt(args[1].substring(8));
                }
                qc.run(mapKey);
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            // parse refseq files from NCBI ftp site
            else if (args[0].contains("refseq")) {

                RefSeqValidator validator = (RefSeqValidator) (bf.getBean("refSeqValidator"));
                validator.run();
                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            // load gene-to-gene associations into RGD_ASSOCIATIONS table
            else if (args[0].contains("load_gene_associations")) {

                int speciesTypeKey = manager.getSpecies(args, 1);

                GeneRelationships relManager = (GeneRelationships) bf.getBean("geneRelationships");
                relManager.loadAssociations(speciesTypeKey);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
            }
            // download all genes for given species
            else if (args[0].contains("all_genes") && args.length>=3) {
                manager.initDbLog(manager.getSpecies(args, 1), "all_genes", args[2]);

                String fileName = manager.entrezGeneExtractor.downloadAndProcessGeneList("all_genes", 
                        EGDAO.getInstance().getEgIdsForAllActiveGenes(manager.speciesTypeKey));
                manager.parseEntrezGeneFile(fileName, null);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
                manager.geneStatusIssueTracker.writeIssuesToFile();
            }
            // download all genes for given species
            else if (args[0].contains("eg_ids_from_file") && args.length>=4) {
                manager.initDbLog(manager.getSpecies(args, 2), "eg_ids_from_file", args[3]);

                String fileName = manager.entrezGeneExtractor.downloadAndProcessGeneList("eg_ids_from_file",
                        EGDAO.getInstance().getEgIdsForActiveGenesFromFile(manager.speciesTypeKey, args[1]));
                manager.parseEntrezGeneFile(fileName, null);

                runSec = (System.currentTimeMillis()-startMilisec)/1000;
                manager.printStatistics();
                manager.logStatistics();
                manager.geneStatusIssueTracker.writeIssuesToFile();
            }
            else {
                manager.printHelp();
                return;
            }              

            manager.logger.info("The process ran "+runSec + " seconds.");
            manager.dbLogger.close(true); // close pipeline log with 'success' indication

            } catch(Exception e) {
                Utils.printStackTrace(e, manager.getLogger());

                // try to write error message and stop db logger
                try {
                    manager.dbLogger.log(e.getMessage(), PipelineLogger.ERROR);
                    manager.dbLogger.close(false); // close pipeline log with 'failure' indication
                }
                catch( Exception ignored ) {
                    ignored.printStackTrace();
                }
                return;
            }
        }
        else {
            manager.printHelp();
            return;
        }   
        manager.getLogger().info("--------------------------------------------------------------------------");
    }
   
    public DataLoadingManager() {
        _instance = this;
    }

    // initialize db logger
    void initDbLog(int species, String mode, String modeInfo) throws Exception {

        dbLogger.init(species, mode, PipelineLogger.PIPELINE_ENTREZGENE);
        dbLogger.log("mode: -"+mode, modeInfo, PipelineLogger.INFO);
        registerDbLogFlags();
    }

    /*
     * First download the data into a file, then chop the file, 
     * convert the file and run the loading
     */
    public void run(String dateFrom, String dateTo) throws Exception {   

        // download external XML file to a local file
        int recordCount = download(dateFrom, dateTo);

        // log info about record count found in the source file
        dbLogger.log("Count of records found in XML file", Integer.toString(recordCount), PipelineLog.LOGPROP_RECCOUNT);

        if (recordCount>0) {
            runEntrezGeneFile(entrezGeneExtractor.getOutFile());            
        }
        else {
            logger.info("There is no data updated during the specified dates"); 
        }
    }
    
    /*
     * Dowload the entrezgene data from NCBI
     * input: from date and outdate (yyyy/mm/dd)
     * output: data file in the specified directory
     * return: nr of records in the local file
     */
    public int download(String dateFrom, String dateTo) throws Exception {

        // handle 'auto' dates
        if( dateFrom.toLowerCase().equals("auto") ) {
            String lastSyncDate = dbLogger.getLastDataSyncDate();
            if( lastSyncDate==null ) {
                // last data sync date not found for pipeline in the current run mode
                // set it to a month ago from the current date
                dateFrom = Utils.addDaysToDate((String)null, -30);
            }
            else
                dateFrom = Utils.addDaysToDate(lastSyncDate, 1);
        }

        String yesterdayDate = Utils.addDaysToDate((String)null, -1);
        if( dateTo.toLowerCase().equals("auto") ) {
            dateTo = yesterdayDate;
            // by default, we pull data up to yesterday
        }

        // handle 'fullload'
        if( dateFrom.toLowerCase().equals("fullload") ) {
            dateFrom = "2000/01/01";
            this.entrezGeneExtractor.setForceFullLoad(true);
        }

        // validity check: dateFrom must be <= dateTo
        if( dateFrom.compareTo(dateTo)>0 ) {
            // swap the dates
            String dt = dateFrom;
            dateFrom = dateTo;
            dateTo = dt;  // now dateFrom <= dateTo
        }

        // dateFrom must be BEFORE today's date
        if( dateFrom.compareTo(yesterdayDate)>0 ) {
            throw new Exception("the starting date "+dateFrom+" must be earlier than today");
        }
        // if dateTo is set to today or after, set it to yesterday
        if( dateTo.compareTo(yesterdayDate)>0 ) {
            dateTo = yesterdayDate;
        }

        String species = SpeciesType.getTaxonomicName(speciesTypeKey);
        entrezGeneExtractor.setSpecies(species);
        dbLogger.log("download: species", species, PipelineLogger.INFO);

        int recordCount;
        entrezGeneExtractor.setDateFrom(dateFrom);
        entrezGeneExtractor.setDateTo(dateTo);
        dbLogger.log(dateFrom, dateTo, PipelineLog.LOGPROP_DATERANGE);

        recordCount = entrezGeneExtractor.run(this.egIds);
        // download finished successfully

        // save dateTo variable as db property LAST-DATA-SYNC-DATE
        dbLogger.setLastDataSyncDate(dateTo);
        dbLogger.writeLogProps();

        return recordCount;
    }
    
    /*
     * convert an entrezgene file into array of BulkGene objects
     * and analyze them
     */
    public void runEntrezGeneFile(String fileName) throws Exception{

        parseEntrezGeneFile(fileName, null);
    }

    /*
     * transform an EntrezGene XML file into BulkGene XML file
     */
    public void transform(String egFileName, String bgFileName) throws Exception {

        parseEntrezGeneFile(egFileName, bgFileName);
    }

    /*
     * Load the data from entrezgene xml file, parse it one record at a time;
     * when the record is complete, put it into a queue for quality checking
     * input: the bulk gene file name
     */
    public void parseEntrezGeneFile(String fileName, String bulkGeneFileName) throws Exception {

        CounterPool counters = new CounterPool();

        String speciesName = SpeciesType.getCommonName(speciesTypeKey);
        logger.info("Processing entrezgene file:" + fileName);
        dbLogger.log("Starting XML processing of entrezgene file", fileName, PipelineLogger.INFO);

        Map<String, String> genomicAssembliesForSpecies = getGenomicAssembliesForCurrentSpecies();
        qualityCheck.setGenomicAssemblies(genomicAssembliesForSpecies);
        qualityCheck.setDbFlagManager(dbFlagManager);
        qualityCheck.setCounters(counters);

        decisionMaker.setDbFlagManager(dbFlagManager);

        // setup thread for concurrent xml parsing
        XomEntrezGeneAnalyzer parser = new XomEntrezGeneAnalyzer();
        parser.setFirstRecNo(getFirstRecNo());
        parser.setGenomicAssemblies(genomicAssembliesForSpecies, getScaffoldAssemblies());
        parser.setGeneLocationHistory(getGeneLocationHistoryForCurrentSpecies());
        parser.setFileName(fileName);

        List<BulkGene> bulkGenes = parser.process(counters);

        // generate file name for storing bulkgenes
        String format = "data/BulkGenes_%s_%s_%s.xml";
        String bgFileName = bulkGeneFileName;
        if( bgFileName==null ) {
            if( entrezGeneExtractor.getDateFrom()==null || entrezGeneExtractor.getDateTo()==null ) {
                String mode = this.dbLogger.getPipelineRunMode();
                SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd_HHmmss");
                String timeNow = sdf.format(new Date());

                bgFileName = String.format(format,
                    speciesName,
                    mode,
                    timeNow);
            }
            else {
                bgFileName = String.format(format,
                    speciesName,
                    this.entrezGeneExtractor.getDateFrom().replace("/",""),
                    this.entrezGeneExtractor.getDateTo().replace("/",""));
            }
        }

        // create one data loading thread
        DataLoadingThread dl = new DataLoadingThread(bgFileName, decisionMaker);
        dl.open();

        for( BulkGene bg: bulkGenes ) {
            qualityCheck.process(bg);
            dl.process(bg, counters);
            counters.increment("GENES_PROCESSED");
        }
        dl.close();

        // store count of records processed
        bgTotal = bulkGenes.size();
        genesWithWrongType = counters.get("WRONG_GENE_TYPE");

        dbLogger.log("Finished XML processing of entrezgene file", fileName, PipelineLogger.INFO);

        // dump counts for genomic assembly names
        for(Map.Entry<String, Integer> entry: parser.getGenomicAssemblyNameCountMap().entrySet() ) {
            // also print out the map key
            //dbLogger.log("[ACTIVE] Assembly Count: "+entry.getKey()+" [MAP_KEY:"+getGenomicAssemblies().get(entry.getKey())+"]", entry.getValue().toString(), PipelineLogger.TOTAL);
            dbLogger.log("[ACTIVE] Assembly Count: "+entry.getKey(), entry.getValue().toString(), PipelineLogger.TOTAL);
        }

        // dump counts for all assemblies found in source data
        Set<String> ignoredAssemblies = new HashSet<String>(parser.getAnyAssemblyNameCountMap().keySet());
        ignoredAssemblies.removeAll(parser.getGenomicAssemblyNameCountMap().keySet());
        for( String ignoredAssembly: ignoredAssemblies ) {
            dbLogger.log("[IGNORED] Assembly Count: "+ignoredAssembly, parser.getAnyAssemblyNameCountMap().get(ignoredAssembly).toString(), PipelineLogger.TOTAL);
            System.out.println("WARNING! Positions found for assembly ["+ignoredAssembly+"] : "+ parser.getAnyAssemblyNameCountMap().get(ignoredAssembly).toString());
        }

        TranscriptVersionManager.getInstance().qcAndLoad(counters);

        // dump counter statistics
        Enumeration<String> counterNames = counters.getCounterNames();
        Set<String> sortedCounterNames = new TreeSet<>();
        while( counterNames.hasMoreElements() ) {
            String counter = counterNames.nextElement();
            sortedCounterNames.add(counter);
        }
        for( String counter: sortedCounterNames ) {
            int count = counters.get(counter);
            if( count!=0 ) {
                count = Math.abs(count);
                dbLogger.log(counter, Integer.toString(count), PipelineLogger.TOTAL);
                System.out.println(counter +": "+count);
            }
        }

        // dump list of eg-ids missed by the pipeline
        String dump = dumpEgIds();
        if( dump!=null ) {

            // dump the eg id list as a string
            dbLogger.getPipelineLogDAO().logDbMessage("missing EG-ID s", null, PipelineLogger.ERROR, dump, dbLogger.getPipelineLog());
        }

        System.out.println("---PROCESSING COMPLETE--");
    }

    private String dumpEgIds() {

        if( egIds.size() > 0 ) {
            StringBuilder buf = new StringBuilder("missing EntrezGene IDs:\n");
            for(Map.Entry<Integer,Integer> entry: egIds.entrySet()) {
                buf.append('[').append(entry.getValue()).append("] ").append(entry.getKey()).append('\n');
            }
            return buf.toString();
        }
        else
            return null;
    }

    /*
     * Print some statistics data for the loading
     */
    public void printStatistics() throws Exception {
        logger.info("Total "+ bgTotal + " entrezgenes are processed.");
        logger.info(genesWithWrongType + " entrezgenes are skipped because of wrong type.");
        if (speciesTypeKey==SpeciesType.HUMAN || speciesTypeKey==SpeciesType.MOUSE) {
            logger.info(decisionMaker.getNoMHID() + " entrezgenes are skipped because they have no MGI/HGNC ID.");
        }
        logger.info(decisionMaker.getNewGenes() + " entrezgenes are new and inserted.");
        logger.info(decisionMaker.getUpdated() + " entrezgenes are updated.");        
        logger.info(decisionMaker.getError() + " entrezgenes have conflicts and needs to be manually checked.");
        logger.info(decisionMaker.skipped + " entrezgenes are skipped (DISCONTINUED, SECONDARY, or gene type is 'other' or 'uknown')");

        // db logging for totals as well
        dbLogger.log("totalNumberGenes: total entrezgene records processed", Integer.toString(bgTotal), PipelineLogger.TOTAL);
        if (speciesTypeKey==SpeciesType.HUMAN || speciesTypeKey==SpeciesType.MOUSE) {
            dbLogger.log("genesNoMHID: entrezgene records skipped because they have no MGI/HGNC ID", Integer.toString(decisionMaker.getNoMHID()), PipelineLogger.TOTAL);
        }
        dbLogger.log("genesNew: entrezgene records new and inserted", Integer.toString(decisionMaker.getNewGenes()), PipelineLogger.TOTAL);
        dbLogger.log("genesUpdated: entrezgene records updated", Integer.toString(decisionMaker.getUpdated()), PipelineLogger.TOTAL);
        dbLogger.log("genesSkipped: (DISCONTINUED, SECONDARY or egType=('other','unknown')", Integer.toString(decisionMaker.skipped), PipelineLogger.TOTAL);
        dbLogger.log("genesNameMismatch: entrezgene records with name different than rgd name", Integer.toString(decisionMaker.getNameMismatch()), PipelineLogger.TOTAL);
        dbLogger.log("timeToExecute: pipeline run time in seconds", Long.toString(runSec), PipelineLogger.TOTAL);
    }
    
    /*
     * Log the loading statistic data into RGD database
     */
    public void logStatistics() throws Exception {

        // log the total records processed
        String subSystem=null;
        if (speciesTypeKey==SpeciesType.RAT) subSystem="entrezGeneRat";
        else if (speciesTypeKey==SpeciesType.MOUSE) subSystem="entrezGeneMouse";
        else if (speciesTypeKey==SpeciesType.HUMAN) subSystem="entrezGeneHuman";

        if (subSystem != null) {            
            rgdLogger.log(subSystem, "totalNumberGenes", bgTotal);
            rgdLogger.log (subSystem, "timeToExecute", Long.toString(runSec));
            rgdLogger.log (subSystem, "genesWrongType", genesWithWrongType);
            rgdLogger.log (subSystem, "genesNew", decisionMaker.getNewGenes());
            rgdLogger.log (subSystem, "genesUpdated", decisionMaker.getUpdated());
            rgdLogger.log (subSystem, "genesSkipped", decisionMaker.skipped);
            if (subSystem.equals("entrezGeneMouse") || subSystem.equals("entrezGeneHuman")) {
                rgdLogger.log(subSystem, "genesNoMHID", decisionMaker.getNoMHID());
            }
        }        
    }
    
    /*
     * Print the information on how to run this pipeline
     */
    public void printHelp() {
        String info =
            "Please enter command line option:\n" +
            " [-download+process | -download_only | -process_only | -refseq | -transform | -load_gene_associations \n"+
            "  | -delete_redundant_aliases] arg1 arg2 ...\n" +
            " -species [rat|mouse|human]\n" +
            "---------------------------------\n" +
            "-download+process dateFrom dateTo\n" +
            "      load data from NCBI(first download the data, then load it)\n" +
            "-download_only outputFileName dateFrom dateTo\n" +
            "      download EntrezGene XML file from NCBI\n" +
            "-process_only entrezGeneFileName\n" +
            "      load data from a specified EntrezGene XML file\n" +
            "-transform inputFileName outputFileName\n" +
            "     transform EntrezGene XML file to BulkGene\n" +
            "-refseq\n" +
            "     generate reports about RefSeq coverage\n" +
            "-refseq_removed\n" +
            "     compare list of removed RefSeq acc ids against RGD and generate a report\n" +
            "-refseq_qc_protein -mapKey=<map-key>\n" +
            "     compare proteins translated from a transcript against RefSeq proteins for given assembly\n" +
            "-load_gene_associations\n" +
            "     download gene_group.gz file from NCBI and load all gene relationships into RGD_ASSOCIATIONS table\n" +
            "-mitochondrial\n" +
            "      download and process all mitochondrial genes from NCBI\n" +
            "-microRNA\n" +
            "      download and process all microRNA genes found in RGD\n" +
            "\n"+
            "     'dateFrom' and 'dateTo' must be formatted as 'yyyy/mm/dd'\n" +
            "     if 'dateFrom' is 'auto', it is set to the 'dateTo' date of latest successful pipeline run\n" +
            "        however, if there was no successful pipeline run yet, it is set to a month ago from today\n" +
            "     if 'dateTo' is 'auto', it is set to yesterday's date\n" +
            "     if 'dateFrom' is 'fullload', no dates are used, and all genes for given species are downloaded and processed\n";
        
        logger.error(info);
        System.out.println(info);

       }

    // process command line arguments and try to enforce new species type key ;
    // return new species type key as read from cmd line or negative value on error
    int getSpecies(String[] args, int argIndex) throws Exception {
        // make sure there is two more arguments available
        if( args.length <= argIndex+1 ) {
            throw new Exception("'-species' argument not found on command line");
        }
        // first argument must be '-species' string literal
        if( !args[argIndex++].equals("-species") ) {
            throw new Exception("'-species' argument not found on command line");
        }

        // 2nd argument must be a valid species type string or key
        int species = SpeciesType.parse(args[argIndex]);
        if( species>0 ) {
            this.speciesTypeKey = species; // species type key override
            BulkGene.setPrimaryMapKey(EGDAO.getInstance().getPrimaryMapKey(species));
            return species;
        }
        else {
            throw new Exception("'-species' argument must be either 'human', 'mouse' or 'rat'");
        }
    }

    // get genomic assemblies to be imported by the pipeline for the currently chosen species
    public Map<String, String> getGenomicAssembliesForCurrentSpecies() {

        String speciesName = SpeciesType.getShortName(speciesTypeKey);
        return getGenomicAssemblies().get(speciesName);
    }

    public DecisionMaker getDecisionMaker() {
        return decisionMaker;
    }
    public void setDecisionMaker(DecisionMaker decisionMaker) {
        this.decisionMaker = decisionMaker;
    }
    public QualityCheckBulkGene getQualityCheck() {
        return qualityCheck;
    }
    public void setQualityCheck(QualityCheckBulkGene qualityCheck) {
        this.qualityCheck = qualityCheck;
    }
    
    public EntrezGeneExtractor getEntrezGeneExtractor() {
        return entrezGeneExtractor;
    }
    public void setEntrezGeneExtractor(EntrezGeneExtractor entrezGeneExtractor) {
        this.entrezGeneExtractor = entrezGeneExtractor;
    }
    
    public RGDSpringLogger getRgdLogger() {
        return rgdLogger;
    }
    public void setRgdLogger(RGDSpringLogger rgdLogger) {
        this.rgdLogger = rgdLogger;
    }
    public Logger getLogger() {
        return this.logger;
    }    
    
    public int getFirstRecNo() {
        return firstRecNo;
    }

    public void setFirstRecNo(int firstRecNo) {
        this.firstRecNo = firstRecNo;
    }

    public int getQualityCheckingThreadCount() {
        return qualityCheckingThreadCount;
    }

    public void setQualityCheckingThreadCount(int qualityCheckingThreadCount) {
        this.qualityCheckingThreadCount = qualityCheckingThreadCount;
    }

    public NcbiEutils geteUtils() {
        return eUtils;
    }

    public void seteUtils(NcbiEutils eUtils) {
        this.eUtils = eUtils;
    }

    public Map<Integer, Integer> getEgIds() {
        return egIds;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    void registerDbLogFlags() throws Exception {
        dbFlagManager.registerFlag(
            "TRANSCRIPT_REFSEQ_STATUS_CHANGED",
            "RefSeq Status for a transcript has been changed"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_PEPTIDE_LABEL_CHANGED",
            "peptide label for a transcript has been changed"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_CODING_STATUS_CHANGED",
            "coding status for a transcript has been changed"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_PROTEIN_ACC_ID_CHANGED",
            "protein acc id for a transcript has been changed"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_INSERTED",
            "new transcript has been inserted"
            );
        dbFlagManager.registerFlag(
            "GENE_MAPPOS_INSERTED",
            "new gene position has been inserted"
            );
        dbFlagManager.registerFlag(
            "GENE_MAPPOS_UPDATED",
            "gene position has been updated"
            );
        dbFlagManager.registerFlag(
            "GENE_MAPPOS_DELETED",
            "gene position has been deleted"
            );
        dbFlagManager.registerFlag(
            "GENE_MAPPOS_DELETE_SUPPRESSED",
            "suppressed deletion of a gene position"
        );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_MAPPOS_INSERTED",
            "new transcript position has been inserted"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_MAPPOS_UPDATED",
            "transcript position has been updated"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_MAPPOS_DELETED",
            "transcript position has been deleted"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_MAPPOS_DELETE_SUPPRESSED",
            "suppressed deletion of transcript position"
            );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_OVERLAPPING_MAPPOS_DELETED",
            "deletion of suppressed transcript position that overlaps other transcript position"
            );
        dbFlagManager.registerFlag(
                "TRANSCRIPT_OVERLAPPING_MAPPOS_DELETE_SUPPRESSED",
                "suppressed deletion of suppressed transcript position that overlaps other transcript position"
        );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_DETACHED_FROM_GENE",
            "transcript detached from gene because it is no longer present in the incoming data"
        );
        dbFlagManager.registerFlag(
            "TRANSCRIPT_DETACH_FROM_GENE_SUPPRESSED",
            "old transcripts have been preserved even if there are no transcripts available for the current assembly"
        );

        dbFlagManager.registerFlag(
            "GENE_TRACK_STATUS_DISCONTINUED",
            "gene record is discontinued - never to be added into rgd"
            );
        dbFlagManager.registerFlag(
            "GENE_TRACK_STATUS_REPLACED",
            "gene record is no longer current because it has been made secondary to another gene record"
            );
        dbFlagManager.registerFlag(
            "GENE_TRACK_STATUS_LIVE",
            "gene record is current and primary"
            );
        dbFlagManager.registerFlag(
            "GENE_TRACK_STATUS_UNKNOWN",
            "gene track status information is not available"
            );
        dbFlagManager.registerFlag(
            "GENE_TRACK_STATUS_SECONDARY",
            "gene record is secondary"
            );

        dbFlagManager.registerFlag(
            "ORTHO_NOMEN_SYMBOL",
            "gene symbol changed (nomenclature event)"
            );
        dbFlagManager.registerFlag(
            "ORTHO_NOMEN_NAME",
            "gene name changed (nomenclature event)"
            );
        dbFlagManager.registerFlag(
            "ORTHO_NOMEN_DESC",
            "gene description changed"
            );
        dbFlagManager.registerFlag(
            "ORTHO_NOMEN_SYMBOL_CHANGE_SUPPRESSED",
            "gene symbol in RGD different from NCBI; change suppressed by HGNC authority"
        );
        dbFlagManager.registerFlag(
            "ORTHO_NOMEN_NAME_CHANGE_SUPPRESSED",
            "gene name in RGD different from NCBI; change suppressed by HGNC authority"
        );

        dbFlagManager.registerFlag(
            "BAD_SPECIES",
            "unexpected species"
            );
        dbFlagManager.registerFlag(
            "EG_IN_RGD",
            "incoming EG is associated with one gene in RGD"
            );
        dbFlagManager.registerFlag(
            "WRONG_GENE_TYPE",
            "incoming gene has wrong gene type: other or unknown"
        );
        dbFlagManager.registerFlag(
            "NEW_GENE_SKIPPED",
            "new incoming gene skipped from loading, because 1) it is discontinued or secondary, "+
                    "or 2) its gene type is 'other'/'unknown' and it does not have any of NM_, NR_, XM_, XR_ or NG_ sequences"
        );
        dbFlagManager.registerFlag(
            "NEW_GENE_INSERTED",
            "new incoming gene inserted into RGD"
        );
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setGeneLocationHistory(Map<String,String> geneLocationHistory) {
        this.geneLocationHistory = geneLocationHistory;
    }

    public Map<String,String> getGeneLocationHistory() {
        return geneLocationHistory;
    }

    public Map<String,Integer> getGeneLocationHistoryForCurrentSpecies() throws Exception {
        Map<String,Integer> historyLocMap = new HashMap<String, Integer>();
        for( Map.Entry<String,String> entry: geneLocationHistory.entrySet() ) {
            int mapKey = Integer.parseInt(entry.getValue());
            if( MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey()==this.speciesTypeKey ) {
                historyLocMap.put(entry.getKey(), mapKey);
            }
        }
        return historyLocMap;
    }

    public GeneStatusIssueTracker getGeneStatusIssueTracker() {
        return geneStatusIssueTracker;
    }

    public void setGenomicAssemblies(Map<String,Map<String,String>> genomicAssemblies) {
        this.genomicAssemblies = genomicAssemblies;
    }

    public Map<String,Map<String,String>> getGenomicAssemblies() {
        return genomicAssemblies;
    }

    public void setScaffoldAssemblies(Map<Integer, String> scaffoldAssemblies) {
        this.scaffoldAssemblies = scaffoldAssemblies;
    }

    public Map<Integer, String> getScaffoldAssemblies() {
        return scaffoldAssemblies;
    }
}
