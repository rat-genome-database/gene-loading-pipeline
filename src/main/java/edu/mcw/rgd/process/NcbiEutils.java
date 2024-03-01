package edu.mcw.rgd.process;

import nu.xom.Builder;
import nu.xom.Document;
import nu.xom.Element;
import org.jaxen.XPath;
import org.jaxen.xom.XOMXPath;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author mtutaj
 * @since 10/18/11
 * A wrapper to automate calls to NCBI eUtils
 */
public class NcbiEutils {

    private String eUtils_db;   // NCBI: database we are connecting to, like 'gene'
    private String eUtils_tool; // name of rgd tool requesting data from NCBI, like 'rgd_eg_pipeline'
    private String eUtils_email; // email contact used to resolve issue between rgd and eUtils
    private String ncbiSearchUrl;
    private String ncbiFetchUrl;

    public ESearchResult createESearchResult() {

        return new ESearchResult();
    }

    /**
     * run eSearch query
     * <p>Note: eSearch result is an xml file; as of Feb 21, 2014, dtd file referenced in result file
     * is no longer available from NCBI eUtils website, therefore to make the parsing work
     * we remove &lt;!DOCTYPE ...&gt; line from result xml file before proceeding
     * @param term term for eSearch query
     * @return ESearchResult object
     * @throws Exception
     */
    public ESearchResult runESearch(String term) throws Exception  {

        // construct eSearch query
        StringBuilder url = new StringBuilder(this.getNcbiSearchUrl());
        url.append("?db=").append(this.geteUtils_db());
        url.append("&tool=").append(this.geteUtils_tool());
        url.append("&email=").append(this.geteUtils_email());
        url.append("&term=").append(term);
        url.append("&retmax=1000000&usehistory=y");

        ESearchResult res = createESearchResult();
        res.queryOriginal = url.toString();
        Document doc = downloadXmlFileFromWeb(res.queryOriginal);
        res.xml = doc.toXML();

        // extract values of elements important to us (parameters to eFetch)
        XPath xpath = new XOMXPath("/eSearchResult/Count");
        res.recordCount = xpath.stringValueOf(doc);
        xpath = new XOMXPath("/eSearchResult/QueryKey");
        res.queryKey = xpath.stringValueOf(doc);
        xpath = new XOMXPath("/eSearchResult/WebEnv");
        res.webEnv = xpath.stringValueOf(doc);
        xpath = new XOMXPath("/eSearchResult/QueryTranslation");
        res.queryTranslation = xpath.stringValueOf(doc);

        // get all entrezgene ids
        xpath = new XOMXPath("/eSearchResult/IdList/Id");
        List<Element> idElList = (List<Element>) xpath.selectNodes(doc);
        res.ids = new ArrayList<Integer>(idElList.size());
        for( Element idEl: idElList ) {
            res.ids.add(Integer.parseInt(idEl.getValue()));
        }

        return res;
    }

    /**
     * based on results from eSearch, we run eFetch, download the file and return a copy of the file;
     * <p>
     * sometimes batch size (value of retMax) is to big for eUtils; if this is the case, we reduce
     * retMax by factor of 5, and try again; eSearch.fetchedCount shows the actual value of fetched records
     * @param eSearch eSearch result
     * @param retMode return mode; if null 'xml' will be used
     * @param retMax maximum nr of records to be returned
     * @param retStart from which record we should fetch the data
     * @return file with data or exception will be thrown
     */
    public File runEFetch(ESearchResult eSearch, String retMode, int retMax, int retStart) throws Exception {

        while( retMax>0 ) {
            int maxRetryCount = 13;
            int retryCount = 0;
            try {
                Thread.sleep(3333); // sleep 3.3 sec

                // construct efetch query
                StringBuilder url = new StringBuilder(this.getNcbiFetchUrl());
                url.append("?db=").append(this.geteUtils_db());
                url.append("&tool=").append(this.geteUtils_tool());
                url.append("&email=").append(this.geteUtils_email());
                url.append("&WebEnv=").append(eSearch.webEnv);
                url.append("&query_key=").append(eSearch.queryKey);
                url.append("&retmode=").append(retMode!=null?retMode:"xml");
                url.append("&retmax=").append(retMax);
                url.append("&retstart=").append(retStart);

                eSearch.fetchQuery = url.toString();
                File f = downloadFile(maxRetryCount, eSearch);

                // the downloaded file must contain <Entrezgene> tag
                if( invalidFileContents(f) )
                    throw new Exception("invalid file contents");

                eSearch.fetchedCount = retMax;
                return f;
            }
            catch(IOException io) {
                // rethrow IO exceptions
                throw io;
            }
            catch(Exception e) {
                // catch FileDownloader exceptions to reduce retMax value by a factor of 5
                if( retMax > 5 ) {
                    retMax /= 5;
                } else if( retMax > 1 ) {
                    retMax = 1;
                } else {
                    if( retryCount<maxRetryCount ) {
                        retryCount++;
                        Thread.sleep(11111*retryCount);
                    } else
                        retMax = 0; // will break the loop
                }
            }
        }
        throw new Exception("permanent error: cannot download data from NCBI");
    }

    private boolean invalidFileContents(File f) throws Exception  {

        // the file must contain <Entrezgene> tag
        BufferedReader reader = new BufferedReader(new FileReader(f));
        String line;
        while( (line=reader.readLine())!=null ) {
            if( line.contains("<Entrezgene>") ) {
                reader.close();
                return false;
            }
        }
        reader.close();
        return true;
    }

    /**
     * based on results from eSearch, we run eFetch, download the file and return a copy of the file
     * @param eSearch eSearch result
     * @param retMode return mode; if null 'xml' will be used
     * @param idList list of eg ids
     * @return file with data
     */
    public File runEFetch(ESearchResult eSearch, String retMode, String idList, int maxRetryCount) throws Exception {

        // construct efetch query
        StringBuilder url = new StringBuilder(this.getNcbiFetchUrl());
        url.append("?db=").append(this.geteUtils_db());
        url.append("&tool=").append(this.geteUtils_tool());
        url.append("&email=").append(this.geteUtils_email());
        url.append("&retmode=").append(retMode!=null?retMode:"xml");
        url.append("&id=").append(idList);

        eSearch.fetchQuery = url.toString();
        return downloadFile(maxRetryCount, eSearch);
    }

    /**
     * based on results from eSearch, we run eFetch, download the file and return a copy of the file
     * @param eSearch eSearch result
     * @param retMode return mode; if null 'xml' will be used
     * @param idList list of eg ids
     * @return file with data
     */
    public File runEFetch(ESearchResult eSearch, String retMode, String idList) throws Exception {

        int maxRetryCount = 8; // default max retry count
        return runEFetch(eSearch, retMode, idList, maxRetryCount);
    }

    public List<Integer> getRelatedPseudoGenes(String egId) throws Exception  {

        ESearchResult res = runESearch("related_functional_gene_"+egId+"[group]");
        return res.ids;
    }

    private Document downloadXmlFileFromWeb(String url) throws Exception {
        // download the file to a tmp file
        File tmpFile = downloadFile(url);

        // remove !DOCTYPE tag
        removeDocType(tmpFile);

        // build xml doc from that file
        Builder builder = new Builder();
        Document doc = null;

        // sometimes we get a 403 response from NCBI server; after a timeout it is usually normal again
        for( int retries=0; retries<10; retries++) {
            try {
                doc = builder.build(tmpFile);

                // sleep up to 1 second after the download so we won't flood NCBI with multiple requests per second
                // what results in getting 403 HTTP errors
                Thread.sleep(new Random().nextInt(1000));

                break;
            }
            catch(java.io.IOException ioe) {
                System.out.println(ioe.toString()+", processing will be retried in few seconds...");

                // sleep 5sec before 1st retry, 10sec before next retry, until 10th retry
                Thread.sleep(5000*(retries+1));
            }
        }

        // cleanup
        tmpFile.delete();
        return doc;
    }

    // this is a hack since https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd file
    // cannot be download from NCBI
    private void removeDocType(File f) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(f));
        File outFile = File.createTempFile("eutils", "tmp");
        BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
        String line;
        while( (line=reader.readLine())!=null ) {
            if( !line.startsWith("<!DOCTYPE") ) {
                writer.write(line);
                writer.newLine();
            }
        }
        reader.close();
        writer.close();

        reader = new BufferedReader(new FileReader(outFile));
        writer = new BufferedWriter(new FileWriter(f));
        while( (line=reader.readLine())!=null ) {
            writer.write(line);
            writer.newLine();
        }
        reader.close();
        writer.close();
    }

    /**
     * download a file given its url, and save it into 'data' subfolder
     * of current program working directory; name of the file will be randomly generated,
     * having prefix 'tmp' and suffix 'xml'
     * @param url url of file to be downloaded
     * @return name of the file created
     * @throws Exception
     */
    public File downloadFile(String url) throws Exception {

        int maxDownloadRetryCount = 8;
        ESearchResult eSearchResult = new ESearchResult();
        eSearchResult.fetchQuery = url;
        return downloadFile(maxDownloadRetryCount, eSearchResult);
    }

    /**
     * download a file given its url, and save it into 'data' subfolder
     * of current program working directory; name of the file will be randomly generated,
     * having prefix 'tmp' and suffix 'xml'
     * @param maxRetryCount maximum number of attempts made to download the file
     * @return name of the file created
     * @throws Exception
     */
    public File downloadFile(int maxRetryCount, ESearchResult eSearchResult) throws Exception {

        // download the file to a tmp file
        FileDownloader downloader = new FileDownloader();
        downloader.setMaxRetryCount(maxRetryCount);
        downloader.setDownloadRetryInterval(20); // set timeout between next download attempt to 20s

        File tmpFile = File.createTempFile("tmp", "xml", new File("data"));
        eSearchResult.localFile = tmpFile.getAbsolutePath();

        downloader.setExternalFile(eSearchResult.fetchQuery);
        downloader.setLocalFile(eSearchResult.localFile);

        downloader.download();

        return tmpFile;
    }

    public String geteUtils_db() {
        return eUtils_db;
    }

    public void seteUtils_db(String eUtils_db) {
        this.eUtils_db = eUtils_db;
    }

    public String geteUtils_tool() {
        return eUtils_tool;
    }

    public void seteUtils_tool(String eUtils_tool) {
        this.eUtils_tool = eUtils_tool;
    }

    public String geteUtils_email() {
        return eUtils_email;
    }

    public void seteUtils_email(String eUtils_email) {
        this.eUtils_email = eUtils_email;
    }

    public String getNcbiFetchUrl() {
        return ncbiFetchUrl;
    }
    public void setNcbiFetchUrl(String ncbiFetchUrl) {
        this.ncbiFetchUrl = ncbiFetchUrl;
    }
    public String getNcbiSearchUrl() {
        return ncbiSearchUrl;
    }
    public void setNcbiSearchUrl(String ncbiSearchUrl) {
        this.ncbiSearchUrl = ncbiSearchUrl;
    }


    /**
     * holds all information relevant in the result from eSearch query
     */
    public class ESearchResult {

        public String recordCount;
        public String queryKey;
        public String webEnv;
        public String queryTranslation; // actual query returned from eUtils
        public String queryOriginal; // original query sent to eUtils
        public String xml; // raw response from eSearch in xml format
        public List<Integer> ids;
        public String fetchQuery; // filled up by NcbiEutils.runEFetch
        public int fetchedCount; // count of records actually fetched by eFetch query
        public int totalFetchedCount;
        public String localFile; // name of local file where contents downloaded via eUtils is stored
    }
}
