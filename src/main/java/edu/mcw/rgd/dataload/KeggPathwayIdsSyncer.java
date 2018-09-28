package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.pipelines.RgdObjectSyncer;
import edu.mcw.rgd.process.Utils;

import java.util.Date;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 11/13/12
 * Time: 11:48 AM
 */
public class KeggPathwayIdsSyncer extends RgdObjectSyncer {

    public EGDAO dao;

    @Override
    protected boolean equalsByUniqueKey(Object o, Object o1) {

        XdbId x1 = (XdbId) o;
        XdbId x2 = (XdbId) o1;
        return Utils.stringsAreEqual(x1.getAccId(), x2.getAccId());
    }

    @Override
    protected boolean equalsByContents(Object o, Object o1) {

        XdbId x1 = (XdbId) o;
        XdbId x2 = (XdbId) o1;
        return Utils.stringsAreEqual(x1.getAccId(), x2.getAccId())
            && Utils.stringsAreEqual(x1.getSrcPipeline(), x2.getSrcPipeline())
            && Utils.stringsAreEqual(x1.getLinkText(), x2.getLinkText());
    }

    @Override
    protected List getDataInRgd(int rgdId) throws Exception {

        return dao.getXdbIdsByRgdId(XdbId.XDB_KEY_KEGGPATHWAY, rgdId);
    }

    @Override
    protected int insertDataIntoRgd(List list) throws Exception {

        return -1;
    }

    @Override
    protected int updateDataInRgd(List list) throws Exception {

        return -1;
    }

    @Override
    protected int deleteDataFromRgd(List list) throws Exception {

        return -1;
    }

    @Override
    protected void copyObjectUniqueKey(Object o, Object o1) {
        XdbId x1To = (XdbId) o;
        XdbId x2From = (XdbId) o1;
        x1To.setKey(x2From.getKey());
        x1To.setCreationDate(x2From.getCreationDate());
        Date sysDate = new Date();
        x1To.setModificationDate(sysDate);
        // do not overwrite incoming link text

        x1To.setNotes("updated by EntrezGene pipeline at "+sysDate);
    }

    @Override
    protected void prepare(Object o, int rgdId, Object userData, int context) {

        XdbId obj = (XdbId) o;
        obj.setRgdId(rgdId);
        obj.setNotes(userData.toString());
    }

    @Override
    public boolean isUpdatable() {
        return true;
    }

    @Override
    public boolean isDeletable() {
        return true;
    }
}
