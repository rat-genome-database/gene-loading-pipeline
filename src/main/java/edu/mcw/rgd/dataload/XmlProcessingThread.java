package edu.mcw.rgd.dataload;

import edu.mcw.rgd.pipelines.RecordPreprocessor;
import edu.mcw.rgd.xml.XomEntrezGeneAnalyzer;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: May 4, 2010
 * Time: 9:35:27 AM
 * reads xml file, breaks into records, parses the records
 * and every parsed record is then put to the queue for further processing
 */
public class XmlProcessingThread extends RecordPreprocessor {

    private final XomEntrezGeneAnalyzer analyzer;

    public XmlProcessingThread(XomEntrezGeneAnalyzer analyzer) {
        this.analyzer = analyzer;
    }

    @Override
    public void process() throws Exception {

        analyzer.setPipelineSession(getSession());

        // use non-validating xml parser
        analyzer.setValidate(false);

        //System.out.println("xml processing thread started!");
        analyzer.parse(new File(analyzer.getFileName()));
        //System.out.println("xml processing thread finished!");
    }
}
