package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.PipelineLogger;

/**
 * @author mtutaj
 * @since May 4, 2010
 */
public class QualityCheckingThread {

    private PipelineLogger dbLogger = PipelineLogger.getInstance();
    private QualityCheckBulkGene qualityCheck;

    public QualityCheckingThread(QualityCheckBulkGene qc) {
        this.qualityCheck = qc;
    }

    public void process(BulkGene bulkGene) throws Exception {

        //System.out.println(Thread.currentThread().getName()+" quality checking thread started!");

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
