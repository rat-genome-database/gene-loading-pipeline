package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.XdbManager;
import org.apache.commons.collections.CollectionUtils;

import java.util.*;

/**
 * @author mtutaj
 * Date: 4/13/2020
 * list of XdbIds present in incoming data for a gene
 */
public class GeneXdbIds {

    private final Set<XdbId> incoming = new HashSet<>();

    public void addIncoming(XdbId id) {
        incoming.add(id);
    }

    public void qc(int rgdId, CounterPool counters) throws Exception {

        EGDAO dao = EGDAO.getInstance();

        XdbId filter = new XdbId();
        filter.setRgdId(rgdId);
        filter.setSrcPipeline(XdbManager.EG_PIPELINE);
        List<XdbId> xdbIdsInRgd = dao.getXdbIds(filter);

        // determine xdb ids for insertion
        Collection<XdbId> forInsert = CollectionUtils.subtract(incoming, xdbIdsInRgd);

        // determine xdb ids for deletion
        Collection<XdbId> forDelete = CollectionUtils.subtract(xdbIdsInRgd, incoming);

        Collection<XdbId> matching = CollectionUtils.intersection(xdbIdsInRgd, incoming);

        qcHgncIds(forInsert, counters, dao);
        qcKeggPathwayIds(forInsert, counters, dao);

        if( !forInsert.isEmpty() ) {
            dao.insertXdbs(forInsert, "GENE");
            counters.add("XDB_IDS_INSERTED", forInsert.size());
        }

        if( !forDelete.isEmpty() ) {
            dao.deleteXdbIds(forDelete, "GENE");
            counters.add("XDB_IDS_DELETED", forDelete.size());
        }

        int matchingXdbIds = matching.size();
        if( matchingXdbIds!=0 ) {
            dao.updateModificationDate(matching);
            counters.add("XDB_IDS_MATCHING", matchingXdbIds);
        }
    }

    /**
     * HGNC IDs that are to be inserted must follow a special rule, that there should not be another active gene
     * that already has this HGNC assigned; any violations must be reported
     * @param forInsert
     */
    void qcHgncIds(Collection<XdbId> forInsert, CounterPool counters, EGDAO dao) throws Exception {

        // currently only report these ids
        Iterator<XdbId> it = forInsert.iterator();
        while( it.hasNext() ) {
            XdbId id = it.next();
            if( id.getXdbKey()==XdbId.XDB_KEY_HGNC ) {
                // insert HGNC ID only if it is not assigned to other gene
                List<RgdId> idsInRgd = dao.getRGDIdsByXdbId(id.getXdbKey(), id.getAccId());
                for( RgdId r: idsInRgd ) {
                    if( r.getObjectKey()==RgdId.OBJECT_KEY_GENES && r.getObjectStatus().equals("ACTIVE") && r.getRgdId()!=id.getRgdId() ) {
                        counters.increment("*** SUPPRESSED INSERTION OF XDB_KEY_HGNC: ACC="+id.getAccId()+" for RGD:"+id.getRgdId());
                        it.remove();
                        break;
                    }
                }
            }
        }
    }

    void qcKeggPathwayIds(Collection<XdbId> forInsert, CounterPool counters, EGDAO dao) throws Exception {

        // special processing rule for Kegg Pathway Ids:
        //
        // currently only report these pathway ids
        Iterator<XdbId> it = forInsert.iterator();
        while( it.hasNext() ) {
            XdbId id = it.next();
            if( id.getXdbKey()==XdbId.XDB_KEY_KEGGPATHWAY ) {
                // insert HGNC ID only if it is not assigned to other gene
                List<RgdId> idsInRgd = dao.getRGDIdsByXdbId(id.getXdbKey(), id.getAccId());
                for( RgdId r: idsInRgd ) {
                    if( r.getObjectKey()==RgdId.OBJECT_KEY_GENES && r.getObjectStatus().equals("ACTIVE") && r.getRgdId()!=id.getRgdId() ) {
                        counters.increment("*** SUPPRESSED INSERTION OF XDB_KEY_KEGGPATHWAY: ACC="+id.getAccId()+" for RGD:"+id.getRgdId());
                        it.remove();
                        break;
                    }
                }
            }
        }
    }
}
