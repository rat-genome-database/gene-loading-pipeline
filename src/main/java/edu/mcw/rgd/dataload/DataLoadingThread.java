package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogger;
import org.apache.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

/**
 * @author mtutaj
 * @since May 5, 2010
 * all database modification code is operated from here
 */
public class DataLoadingThread {
    private PipelineLogger dbLogger = PipelineLogger.getInstance();
    private String bgFileName;
    private DecisionMaker decisionMaker;
    private final Logger logger = Logger.getLogger("process");
    private BufferedWriter out;

    public DataLoadingThread(String bgFileName, DecisionMaker dm) {
        this.bgFileName = bgFileName;
        this.decisionMaker = dm;
    }

    /**
     * open the file for storing bulk genes
     * @throws IOException
     */
    public void open() throws IOException {
        // open the file for storing bulk genes
        this.out = new BufferedWriter(new FileWriter(bgFileName));
        logger.info("opening "+new File(bgFileName).getAbsolutePath());
        //System.out.println("data loading thread started!");
        // and write xml header
        out.write("<?xml version=\"1.0\" ?>\n<EntrezGeneSet version=\"1.0.1\">\n");
    }

    /**
     * close the bulk gene file and write a status message to db pipeline log
     * @throws Exception
     */
    public void close() throws Exception {
        // terminate the output file
        out.write("</EntrezGeneSet>");
        out.close();

        logger.info("closing "+new File(bgFileName).getAbsolutePath());

        // we write bulk genes once again to the file, possibly with modified information acquired during the run-time
        dbLogger.log("Bulk genes written to file", bgFileName, PipelineLogger.INFO);

        //System.out.println("data loading thread finished!");
    }

    /**
     * process a single BulkGene record
     * @throws Exception
     */
    public void process(BulkGene bulkGene, CounterPool counters) throws Exception {

        //System.out.println("LOAD "+bulkGene.getRecNo()+", qsize:"+threadResult.bgDataLoadingQueue.size());

        // decision making
        decisionMaker.decide(bulkGene, counters);

        // additional logging of matching rgd_id that could be known after decision making is finished
        if( bulkGene.getCustomFlags().getRgdId()>0 ) {
            dbLogger.addLogProp("RGD_ID matching the source record", Integer.toString(bulkGene.getCustomFlags().getRgdId()), bulkGene.getRecNo(), PipelineLogger.REC_MATCHINGRGDID);
            dbLogger.writeLogProps(bulkGene.getRecNo());
            dbLogger.removeAllLogProps(bulkGene.getRecNo());
        }

        // write xml representation to the file
        out.write(bulkGene.getXmlString());
        out.write("\n");

        // remove eg id from 'missing EG IDs map'
        // map of EntrezGene IDS returned by eSearch to 0-based position nr (assigned sequentially during creation);
        // when gene record with given EG ID is successfully processed, its EG ID is removed from map;
        // any left-over EG-IDS mean error either in EG downloading module or in eUtils inconsistency:
        // (eFetch return less records than eSearch)
        Map<Integer,Integer> egIds = DataLoadingManager.getInstance().getEgIds();
        egIds.remove(Integer.parseInt(bulkGene.getEgId()));
    }
}
