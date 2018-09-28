package edu.mcw.rgd.dataload;

import edu.mcw.rgd.pipelines.PipelineRecord;
import edu.mcw.rgd.pipelines.RecordProcessor;
import edu.mcw.rgd.process.PipelineLogger;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: May 4, 2010
 * Time: 9:55:09 AM
 */
public class QualityCheckingThread extends RecordProcessor {

    private PipelineLogger dbLogger = PipelineLogger.getInstance();
    private QualityCheckBulkGene qualityCheck;

    public QualityCheckingThread(QualityCheckBulkGene qc) {
        this.qualityCheck = qc;
    }

    public void process(PipelineRecord pipelineRecord) throws Exception {

        BulkGene bulkGene = (BulkGene) pipelineRecord;

        //System.out.println(Thread.currentThread().getName()+" quality checking thread started!");

        bulkGene.setSession(getSession());

        //System.out.println(Thread.currentThread().getName()+" "+bulkGene.getRecNo()+", qsize:"+inputQueue.size());
        // try to perform quality check
        Flags flag = qualityCheck.check(bulkGene);
        bulkGene.setCustomFlags(flag);

        // 2nd pass
        qualityCheck.afterCheck(bulkGene);

        // record analyzed: write all log props to database
        dbLogger.writeLogProps(bulkGene.getRecNo());
        dbLogger.removeAllLogProps(bulkGene.getRecNo());

        //System.out.println(Thread.currentThread().getName()+" quality checking thread finished!");
    }
}
