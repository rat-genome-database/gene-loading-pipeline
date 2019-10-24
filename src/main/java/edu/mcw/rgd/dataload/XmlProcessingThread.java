package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.xml.XomEntrezGeneAnalyzer;

import java.io.File;
import java.util.List;

/**
 * @author mtutaj
 * @since May 4, 2010
 * reads xml file, breaks into records, parses the records
 * and every parsed record is then put to the queue for further processing
 */
public class XmlProcessingThread {

    private final XomEntrezGeneAnalyzer analyzer;

    public XmlProcessingThread(XomEntrezGeneAnalyzer analyzer) {
        this.analyzer = analyzer;
    }

    public List<BulkGene> process(CounterPool counters) throws Exception {

        analyzer.setCounters(counters);

        // use non-validating xml parser
        analyzer.setValidate(false);

        //System.out.println("xml processing thread started!");
        analyzer.parse(new File(analyzer.getFileName()));
        //System.out.println("xml processing thread finished!");

        return analyzer.bulkGenes;
    }
}
