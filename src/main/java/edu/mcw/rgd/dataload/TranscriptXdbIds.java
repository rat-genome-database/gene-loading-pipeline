package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import org.apache.commons.collections.CollectionUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 9/30/15
 * Time: 3:27 PM
 * <p>
 * handle xdb ids for transcripts; currently only UniProt ids
 */
public class TranscriptXdbIds {

    List<TrXdbId> trXdbIds = new ArrayList<>();

    Collection<XdbId> xdbIdsMatching;
    Collection<XdbId> xdbIdsForInsert;
    Collection<XdbId> xdbIdsForDelete;

    public void addUniProtId(String uniProtId, String uniProtType, String productId) {
        //System.out.println(uniProtId+" "+uniProtType+" "+productId);

        // check if no duplicates
        TrXdbId trXdbId = new TrXdbId(uniProtId, uniProtType, productId);
        for( TrXdbId o: trXdbIds ) {
            if( trXdbId.equalsTo(o) ) {
                return; // duplicate found
            }
        }
        trXdbIds.add(trXdbId);
    }

    public void qc(Transcript tr, BulkGene bg) throws Exception {

        // load uniprot ids from RGD for this transcript
        XdbId filter = new XdbId();
        filter.setSrcPipeline("ENTREZGENE");
        filter.setRgdId(tr.getRgdId());
        filter.setXdbKey(XdbId.XDB_KEY_UNIPROT);
        List<XdbId> xdbIdsInRgd = bg.dao.getXdbIds(filter);

        // find the incoming xdb ids for this transcript
        List<XdbId> xdbIdsIncoming = new ArrayList<>();
        for( TrXdbId trXdbId: trXdbIds ) {
            if( Utils.stringsAreEqual(trXdbId.productId, tr.getProteinAccId() ) ) {
                XdbId xdbId = new XdbId();
                xdbId.setXdbKey(XdbId.XDB_KEY_UNIPROT);
                xdbId.setRgdId(tr.getRgdId());
                xdbId.setAccId(trXdbId.uniProtId);
                xdbId.setSrcPipeline("ENTREZGENE");
                xdbId.setNotes(trXdbId.uniProtType);
                xdbIdsIncoming.add(xdbId);
            }
        }

        xdbIdsMatching = CollectionUtils.intersection(xdbIdsIncoming, xdbIdsInRgd);
        xdbIdsForInsert = CollectionUtils.subtract(xdbIdsIncoming, xdbIdsInRgd);
        xdbIdsForDelete = CollectionUtils.subtract(xdbIdsInRgd, xdbIdsIncoming);
    }

    public void load(BulkGene bg, CounterPool counters) throws Exception {

        if( !xdbIdsMatching.isEmpty() ) {
            bg.dao.updateModificationDate(xdbIdsMatching);
            counters.add("TRANSCRIPT_UNIPROT_IDS_MATCHING", xdbIdsMatching.size());
        }

        if( !xdbIdsForInsert.isEmpty() ) {
            bg.dao.insertXdbs(xdbIdsForInsert, "TRANSCRIPT");
            counters.add("TRANSCRIPT_UNIPROT_IDS_INSERTED", xdbIdsForInsert.size());
        }

        if( !xdbIdsForDelete.isEmpty() ) {
            bg.dao.deleteXdbIds(xdbIdsForDelete, "TRANSCRIPT");
            counters.add("TRANSCRIPT_UNIPROT_IDS_DELETED", xdbIdsForDelete.size());
        }
    }

    class TrXdbId {
        String uniProtId;  // f.e. 'P62630'
        String uniProtType; // 'UniProtKB/Swiss-Prot' or 'UniProtKB/TRembl'
        String productId; // f.e. 'NP_787032'

        public TrXdbId(String uniProtId, String uniProtType, String productId) {
            this.uniProtId = uniProtId;
            this.uniProtType = uniProtType;
            this.productId = productId;
        }

        public boolean equalsTo(TrXdbId o) {
            return Utils.stringsAreEqual(o.uniProtId, this.uniProtId) &&
                    Utils.stringsAreEqual(o.uniProtType, this.uniProtType) &&
                    Utils.stringsAreEqual(o.productId, this.productId);
        }
    }
}
