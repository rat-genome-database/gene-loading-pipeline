package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.*;
import edu.mcw.rgd.xml.XomAnalyzer;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * @author dli, mt
 *
 * @since Sep 13, 2006
 * download gene data in xml format and store it in local file
 * <p>
 * Sep 12, 2012 - changed uTils query to use [Modification+Date] instead of [MODDATE]
 *                which has been discontinued at NCBI
 */

public class EntrezGeneExtractor {
    
    String dateFrom; // = "2006/07/23";
    String dateTo; // = "2006/07/30";    
    String outDir;
    String outFile;// the output file
    String species;
    int maxRecCountInSingleFetch;
    int downloadMaxRetries;
    int delayBetweenFetches;
    private boolean useHistory;
    private boolean forceFullLoad;
    private boolean mitochondrialGenesLoad; // set this to true to load mitochondrial genes

    NcbiEutils.ESearchResult eSearchOverride;

    protected final Logger logger = LogManager.getLogger("process");
    PipelineLogger dbLogger = PipelineLogger.getInstance();

    // download the data into a local disk file; return nr of records to be found in the file
    public int run(Map<Integer,Integer> egIds) throws Exception {

        // eUtils helper class
        NcbiEutils eUtils = DataLoadingManager.getInstance().geteUtils();
        NcbiEutils.ESearchResult eSearch = null;
        String speciesName = SpeciesType.getCommonName(SpeciesType.parse(species));

        SimpleDateFormat sdt = new SimpleDateFormat("yyyyMMdd");
        String todayDate = sdt.format(new Date());

        if( eSearchOverride!=null ) {
            eSearch = eSearchOverride;
        } else {

            // prepare eSearch query
            species = species.replaceAll(" ", "+");
            String query = "(%22"+species+"%22[Organism])";

            if( getMitochondrialGenesLoad() ) {
                query += "+AND+MT[Chromosome]";
            } else {

                // as of Feb 2024, we decided to process all genes, including inactive ones
                //query += "+AND+alive[prop]";

                if( !getForceFullLoad() ) {
                    query += "+AND+(%22"+dateFrom+"%22[Modification+Date]+:+%22"+dateTo+"%22[Modification+Date])";
                }
            }

            // Execute query and parse results for webEnv amd queryKey to be used in later operations
            eSearch = eUtils.runESearch(query);
            String searchAction = eSearch.queryOriginal;
            logger.info("Search Action:"+ searchAction);
            dbLogger.log("NCBI eSearch uri", searchAction, PipelineLogger.INFO);

            // dump the document into separate file
            makeDir(outDir);

            String eSearchFile = outDir + "/" + "EntrezGene_"+speciesName;
            if( getForceFullLoad() ) {
                eSearchFile = outDir + "/fullload_"+speciesName+"_"+todayDate;
            } else
                eSearchFile += "_"+getDateFrom().replaceAll("/", "") +"_"+getDateTo().replaceAll("/", "");
            eSearchFile += "_eSearch.xml";

            FileWriter fw = new FileWriter(eSearchFile);
            fw.append(eSearch.xml);
            fw.close();

            dbLogger.log("NCBI eFetch Count", eSearch.recordCount, PipelineLogger.DEBUG);
            dbLogger.log("NCBI eFetch QueryKey", eSearch.queryKey, PipelineLogger.DEBUG);
            dbLogger.log("NCBI eFetch WebEnv", eSearch.webEnv, PipelineLogger.DEBUG);
            dbLogger.log("NCBI eFetch QueryTranslation", eSearch.queryTranslation, PipelineLogger.DEBUG);
        }


        // check if history data is available
        boolean isHistoryAvailable = !Utils.isStringEmpty(eSearch.queryKey) && !Utils.isStringEmpty(eSearch.webEnv);

        // put all entrezgene ids returned by eSearch to a shared map
        for( Integer idEl: eSearch.ids ) {
            egIds.put(idEl, egIds.size());
        }

        logger.info("Total records: "+ eSearch.recordCount);
        logger.debug("queryKey: "+ eSearch.queryKey);
        logger.debug("webEnv: "+ eSearch.webEnv);
        // if there is no data fetched
        int recCount = Integer.valueOf(eSearch.recordCount);
        if ( recCount == 0) {
            logger.info("There is no data changed between the specified dates");
            return 0;
        }

        // all of the above strings must be non-null and non-empty
        if( Utils.isStringEmpty(eSearch.recordCount) || Utils.isStringEmpty(eSearch.queryTranslation) || (!isHistoryAvailable && eSearch.ids.isEmpty())) {
            // no data received -- error
            throw new Exception("failed to get a document from NCBI eSearch");
        }


        // initialize the file for receiving all of the data
        String outFile = outDir + "/" + "EntrezGene_"+speciesName;
        if( getForceFullLoad() )
            outFile += "_fullload_"+todayDate;
        else if( eSearchOverride!=null )
            outFile += eSearchOverride.queryTranslation;
        else
            outFile += "_"+dateFrom.replaceAll("/", "") +"_"+dateTo.replaceAll("/", "");
        outFile += ".xml.gz";
        setOutFile(outFile);

        initEgWriter();

        // override isHistoryAvailable
        if( !getUseHistory() )
            isHistoryAvailable = false; // turn off history mode
        // note: there is no need to turn on history mode because it will be used if available by default

        // there are limits in NCBI eUtils -- as of April 2010, you cannot download more than 10000 gene records at a time
        if( isHistoryAvailable ) {
            // the variant with history available should be used whenever possible
            // because of better performance
            int retstart = 0;
            int fetchBatchSize = maxRecCountInSingleFetch;
            for( int recCountLeft = recCount; recCountLeft>0;  ) {
                int recCountToFetch = fetchBatchSize;
                if( recCountToFetch > recCountLeft )
                    recCountToFetch = recCountLeft;

                // Fetch the results of the previous query as XML records
                File tmpFile = doEFetch(eUtils, eSearch, recCountToFetch, retstart);
                logger.debug("Fetch Action:"+ eSearch.fetchQuery);
                dbLogger.log("NCBI eFetch uri", eSearch.fetchQuery, PipelineLogger.INFO);

                // sometimes ncbi returns less than we would like to
                if( eSearch.fetchedCount < recCountToFetch ) {
                    // ncbi fetched return much less than requested; assuming this is permanent,
                    // we reduce our fetch batch size to average between requested and actual batch size
                    fetchBatchSize = (eSearch.fetchedCount + recCountToFetch) / 2;
                    if( fetchBatchSize <= 0 )
                        fetchBatchSize = 1; // ensure fetch batch size is at least 1
                }
                else {
                    // ncbi returned exactly amount of records as requested
                    // if this is below our default batch size, bump batch size by 10% + 5
                    fetchBatchSize = 5 + (fetchBatchSize*110)/100;
                    if( fetchBatchSize > maxRecCountInSingleFetch )
                        fetchBatchSize = maxRecCountInSingleFetch; // ensure batch size does not exceed default batch size
                }
                recCountToFetch = eSearch.fetchedCount;

                // read the temporary file and append it to our main eg file
                appendToEgFile(tmpFile);
                // delete the temporary file
                tmpFile.delete();


                recCountLeft -= recCountToFetch;
                retstart += recCountToFetch;

                // sleep between consecutive fetches to avoid overload of NCBI server (try to be nice to them :-) )
                sleepBetweenFetchRequests();
            }
        }
        else if( !eSearch.recordCount.equals("0") ) {

            final int idListSize = getMitochondrialGenesLoad() ? 1 : 500;
            List<Integer> egIdSet = new ArrayList<Integer>(egIds.keySet()); // make a copy of eg id set
            Collections.shuffle(egIdSet);
            eSearch.totalFetchedCount = 0;
            downloadWithoutHistory(egIdSet, eUtils, eSearch, idListSize);
        }
        closeEgWriter();

        logger.info("The data is downloaded into file "+ outFile);
        dbLogger.log("Data from NCBI eFetch saved to file", outFile, PipelineLogger.INFO);
        return recCount;
    }

    public String downloadAndProcessGeneList(String mode, Collection<Integer> egIds) throws Exception {

        int speciesTypeKey = DataLoadingManager.getInstance().getSpeciesTypeKey();
        SimpleDateFormat sdt = new SimpleDateFormat("yyyyMMdd");
        String todayDate = sdt.format(new Date());
        setOutFile("data/"+todayDate+"_"+mode+"_"+SpeciesType.getCommonName(speciesTypeKey).toLowerCase()+".xml.gz");
        initEgWriter();

        // download and process the data
        final int idListSize = 300;
        List<Integer> egIdSet = new ArrayList<>(egIds); // make a copy of eg id set
        Collections.shuffle(egIdSet);
        NcbiEutils eUtils = DataLoadingManager.getInstance().geteUtils();
        NcbiEutils.ESearchResult eSearch = eUtils.createESearchResult();
        eSearch.totalFetchedCount = 0;
        eSearch.recordCount = Integer.toString(egIdSet.size());

        downloadWithoutHistory(egIdSet, eUtils, eSearch, idListSize);
        closeEgWriter();

        return getOutFile();
    }

    void downloadWithoutHistory(List<Integer> egIds, NcbiEutils eUtils, NcbiEutils.ESearchResult eSearch, int idListLen) throws Exception {

        int originalIdListLen = idListLen;

        // no history data available -- all eFetch requests post the GeneID-s explicitly
        //
        while( !egIds.isEmpty() ) {

            // create a string of ids; string length must not exceed the given limit
            Iterator<Integer> it = egIds.iterator();
            StringBuilder idListStringBuilder = new StringBuilder(idListLen+32);
            idListStringBuilder.append(it.next().toString()); // at least one egId must be in the list of egIds for processing
            int processedEgIds = 1;
            it.remove();
            while( it.hasNext() && idListStringBuilder.length()<idListLen ) {
                idListStringBuilder.append(',').append(it.next().toString());
                it.remove();
                processedEgIds++;
            }
            String idList = idListStringBuilder.toString();

            // Fetch the results of the previous query as XML records
            int u;
            File tmpFile;
            int retryCount = idListLen > 100 ? 1 : idListLen > 10 ? 3 : 5;
            for( u=0; u<retryCount; u++ ) {

                int fetchDurationInSeconds;
                try {
                    logger.debug("init eFetch for: " + idList);
                    long time0 = System.currentTimeMillis();
                    tmpFile = eUtils.runEFetch(eSearch, "xml", idList, 1);
                    fetchDurationInSeconds = (int) ((System.currentTimeMillis()-time0)/1000);
                }
                catch( Exception e ) {
                    if( e.getMessage().contains("permanent error") ) {
                        // it must be thrown by FileDownloader -- maximum number of download retrials have been reached
                        // try shorter id list
                        u = retryCount + 10;
                        break;
                    }
                    else {
                        // rethrow any other exception
                        throw e;
                    }
                }

                if( verifyTmpFile(tmpFile) ) {
                    // read the temporary file and append it to our main eg file
                    appendToEgFile(tmpFile);

                    // delete the temporary file
                    tmpFile.delete();

                    eSearch.totalFetchedCount += processedEgIds;

                    logger.debug("processed " + processedEgIds + " ids, total processed " + eSearch.totalFetchedCount + " ids out of " + eSearch.recordCount
                            + "\nId List Len : " + idListLen
                            + "\nFetch Action:" + eSearch.fetchQuery);
                    dbLogger.log("NCBI eFetch uri", eSearch.fetchQuery, PipelineLogger.INFO);

                    // sleep between consecutive fetches to avoid overload of NCBI server (try to be nice to them :-) )
                    sleepBetweenFetchRequests();

                    // successful download -- if time duration needed to download the file via eFetch is more than 1 min
                    //    then decrease the list len by 10%
                    if( fetchDurationInSeconds>60 )
                        idListLen -= 1+(idListLen/20);

                    // successful download -- bump up slightly the idListLen
                    if( idListLen < originalIdListLen )
                        // bump it up by 5%
                        idListLen += 1 + (idListLen/20);
                    else
                        idListLen += 3;
                    logger.debug("New Id List Len : " + idListLen);
                    break;
                }
                else {
                    // file verification failed -- sleep and try again download the file
                    logger.warn("download: file verification failed - refetching ..."+
                        "\nId List Len : "+ idListLen);

                    // sleep between consecutive fetches to avoid overload of NCBI server (try to be nice to them :-) )
                    sleepBetweenFetchRequests();
                }
            }

            if( u >= retryCount ) {
                Set<Integer> egIdsProcessedPartially = appendPartialFile(eSearch.localFile);

                // add again all of the ids to the list for download
                int unprocessedEgIds = 0;
                for( String egId: idList.split("[,]") ) {
                    int idToLoad = Integer.parseInt(egId);
                    if( !egIdsProcessedPartially.contains(idToLoad) ) {
                        egIds.add(idToLoad);
                        unprocessedEgIds++;
                    }
                }
                logger.warn("-- "+unprocessedEgIds+" unprocessed Eg Ids returned to id-list");

                processedEgIds -= unprocessedEgIds;
                eSearch.totalFetchedCount += processedEgIds;

                logger.info("-- processed " + processedEgIds + " ids, total processed " + eSearch.totalFetchedCount + " ids out of " + eSearch.recordCount);

                idListLen /= 2;
                if( idListLen <= 0 )
                    idListLen = 1;

                // all download attempts failed -- reduce count of ids on the cmd line and retry download
                logger.warn("download: all attempts failed - retrying with shorter id-list..."+
                        "\nNew Id List Len : "+ idListLen);
                Collections.shuffle(egIds); // randomize id list to be refetched
            }
        }
    }

    // try to download the file through eFetch, 3 times
    private File doEFetch( NcbiEutils eUtils, NcbiEutils.ESearchResult eSearch, int recCountToFetch, int retstart) throws Exception {

        // try 3 times; retrying if IOException happens
        for( int i=0; i<3; i++ ) {
            try {
                return eUtils.runEFetch(eSearch, "xml", recCountToFetch, retstart);
            }
            catch(IOException io) {
                io.printStackTrace(System.out);
            }
        }

        // no mercy now: any exception will terminate the processing
        return eUtils.runEFetch(eSearch, "xml", recCountToFetch, retstart);
    }

    private Writer _egWriter = null;

    // initialize the output file with EntrezGene records and write the header
    private void initEgWriter() throws IOException {

        // open the writer overwriting any previous contents
        _egWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile))));

        // write the header to the file
        String header = "<?xml version=\"1.0\"?>\n" +
                "<Entrezgene-Set>\n\n";
        _egWriter.write(header);
    }

    // the file read may not contain
    private boolean verifyTmpFile (File reader) throws IOException {

        // read the temporary file and append it to our main eg file
        XomAnalyzer analyser = new XomAnalyzer();
        try {
            analyser.parse(reader);
        }
        catch( Exception e) {
            // if we got this exception, the file could not be parsed at all
            return false;
        }
        return true;
    }

    // write to EntrezGene file the contents read from the reader, except the EG header
    private void appendToEgFile (File reader ) throws IOException {

        // read the temporary file and append it to our main eg file
        BufferedReader lineReader = new BufferedReader(new FileReader(reader));

        // read the lines until you find the line with <Entrezgene> element
        String line;
        while( (line=lineReader.readLine())!=null ) {
            if(line.contains("<Entrezgene>")) {
                int pos = line.indexOf("<Entrezgene>");
                _egWriter.write(line.substring(pos));
                _egWriter.append('\n');
                break;
            }
        }

        // read the rest of lines; but if the line contains terminating <Entrezgene-Set>, do NOT write it
        while( (line=lineReader.readLine())!=null ) {
            if(line.contains("</Entrezgene-Set>")) {
                break;
            }
            _egWriter.write(line);
            _egWriter.append('\n');
        }

        // close the reader
        lineReader.close();
    }

    // write the terminating tag to the file and close the file
    private void closeEgWriter() throws IOException {

        _egWriter.write("</Entrezgene-Set>\n");
        _egWriter.close();
        _egWriter = null;
    }

    // the file most likely is a partial (broken) xml file
    // try to extract the whole EntrezGene records from the file
    private Set<Integer> appendPartialFile(String partialFile) throws Exception {
        Set<Integer> egIds = new TreeSet<Integer>();

        String partialFileContents = Utils.readFileAsString(partialFile);

        File outFile = File.createTempFile("partial", "xml", new File("data"));
        FileWriter out = new FileWriter(outFile);

        final String tagStart = "<Entrezgene>";
        final String tagEnd = "</Entrezgene>";
        int offset = 0;
        while( true ) {
            int pos1 = partialFileContents.indexOf(tagStart, offset);
            int pos2 = partialFileContents.indexOf(tagEnd, offset);
            if( pos1<0 || pos2<0 || pos2<=pos1 )
                break;
            pos2 += tagEnd.length();

            // extract whole <Entrezgene> object from partial file
            String tagBody = partialFileContents.substring(pos1, pos2);

            // extract eg id from tag body
            int pos3 = tagBody.indexOf("<Gene-track_geneid>");
            int pos4 = tagBody.indexOf("</Gene-track_geneid>");
            if( pos3>0 && pos4>pos3 ) {
                String egId = tagBody.substring(pos3+19, pos4);
                egIds.add(Integer.parseInt(egId));
                out.append(tagBody);
                offset = pos2;
            } else {
                logger.warn("--PARTIAL FILE problem -- no GeneId");
                break;
            }
        }

        logger.warn("--PARTIAL FILE extracted genes: "+Utils.concatenate(egIds, ","));

        //delete contents of original partial file
        new File(partialFile).delete();

        out.close();

        appendToEgFile(outFile);

        // delete well-formed file with fully extracted genes
        outFile.delete();

        return egIds;
    }

    private boolean makeDir(String directoryName) {
        // Create a directory; all non-existent ancestor directories are
       // automatically created
       boolean success=false;
       boolean exists = (new File(directoryName)).exists();
       if (exists) return success;
       try {
           success = (new File(directoryName)).mkdirs();
           if (!success) {
               logger.error("Making directory "+ directoryName +" failed");
           } 
       } catch (java.lang.SecurityException se) {
           logger.error("Making directory SecurityException catched");
       }
       return success;
    }

    void sleepBetweenFetchRequests() throws InterruptedException {

        // during work hours: 8am through 6 pm, and Monday through Friday, double the sleep duration
        long millis = 1000*getDelayBetweenFetches();

        Calendar cal = new GregorianCalendar();
        int dayOfWeek = cal.get(Calendar.DAY_OF_WEEK);
        int hourOfDay = cal.get(Calendar.HOUR_OF_DAY);
        if( dayOfWeek!=Calendar.SATURDAY && dayOfWeek!=Calendar.SUNDAY && hourOfDay>=8 && hourOfDay<=18 ) {
            millis *= 2;
        }

        Thread.sleep(millis);
    }

    public String getDateFrom() {
        return dateFrom;
    }
    public void setDateFrom(String dateFrom) {
        this.dateFrom = dateFrom;
    }
    public String getDateTo() {
        return dateTo;
    }
    public void setDateTo(String dateTo) {
        this.dateTo = dateTo;
    }
    
    public String getOutFile() {
        return outFile;
    }
    public void setOutFile(String outFile) {
        this.outFile = outFile;
    }
    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }
    
    public String getSpecies() {
        return species;
    }
    public void setSpecies(String species) {
        this.species = species;
    }

    public int getMaxRecCountInSingleFetch() {
        return maxRecCountInSingleFetch;
    }

    public void setMaxRecCountInSingleFetch(int maxRecCountInSingleFetch) {
        this.maxRecCountInSingleFetch = maxRecCountInSingleFetch;
    }

    public int getDownloadMaxRetries() {
        return downloadMaxRetries;
    }

    public void setDownloadMaxRetries(int downloadMaxRetries) {
        this.downloadMaxRetries = downloadMaxRetries;
    }

    public int getDelayBetweenFetches() {
        return delayBetweenFetches;
    }

    public void setDelayBetweenFetches(int delayBetweenFetches) {
        this.delayBetweenFetches = delayBetweenFetches;
    }

    public void setUseHistory(boolean useHistory) {
        this.useHistory = useHistory;
    }

    public boolean getUseHistory() {
        return useHistory;
    }

    public void setForceFullLoad(boolean forceFullLoad) {
        this.forceFullLoad = forceFullLoad;
    }

    public boolean getForceFullLoad() {
        return forceFullLoad;
    }

    public boolean getMitochondrialGenesLoad() {
        return mitochondrialGenesLoad;
    }

    public void setMitochondrialGenesLoad(boolean mitochondrialGenesLoad) {
        this.mitochondrialGenesLoad = mitochondrialGenesLoad;
    }
}

