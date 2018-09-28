package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.NomenclatureEvent;
import edu.mcw.rgd.datamodel.RgdId;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.Timestamp;
import java.util.Date;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 6/29/12
 * Time: 7:57 AM
 */
public class GeneReactivator {

    private String file;
    private String refKey;
    private String nomenEventType;
    private String nomenEventDesc;

    public void run() throws Exception {

        BufferedReader reader = new BufferedReader(new FileReader(getFile()));
        int genesReactivated = 0;
        int genesSkipped = 0;
        String line;
        while( (line=reader.readLine())!=null ) {
            // [2012-06-28 16:29:36,578] - RGDID=5045246, EGID=292228
            // extract RGDID from the line
            int pos1 = line.indexOf('=');
            if( pos1<0 )
                continue;
            int pos2 = line.indexOf(',', pos1);
            if( pos2<0 )
                continue;
            String rgdId = line.substring(pos1+1, pos2);
            if( reactivateGene(Integer.parseInt(rgdId)) ) {
                genesReactivated++;
            }
            else {
                genesSkipped++;
            }
        }
        reader.close();

        System.out.println("Genes reactivated: "+genesReactivated);
        if( genesSkipped>0 ) {
            System.out.println("Genes skipped: "+genesSkipped);
        }
    }

    boolean reactivateGene(int rgdId) throws Exception {

        EGDAO dao = EGDAO.getInstance();

        // make sure that given object is a gene and that it is not active
        RgdId rgd = dao.getRgdId(rgdId);
        if( rgd==null ) {
            System.out.println("Error: "+rgdId+" is not a valid RGD ID!");
            return false;
        }
        if( rgd.getObjectKey()!=RgdId.OBJECT_KEY_GENES ) {
            System.out.println("Error: "+rgdId+" is not a gene! It is "+rgd.getObjectTypeName());
            return false;
        }
        if( rgd.getObjectStatus().equals("ACTIVE") ) {
            System.out.println("Error: "+rgdId+" is active!");
            return false;
        }

        // debug message
        Gene gene = dao.getGene(rgdId);
        System.out.println("Reactivated ["+gene.getType()+"] gene ["+gene.getSymbol()+"] RGDID:"+rgdId+", old status "+rgd.getObjectStatus());

        // create nomenclature event
        NomenclatureEvent event = new NomenclatureEvent();
        event.setRgdId(rgdId);
        event.setRefKey(getRefKey());
        event.setEventDate(new Timestamp(System.currentTimeMillis()));
        event.setSymbol(gene.getSymbol());
        event.setName(gene.getName());
        event.setNomenStatusType(getNomenEventType());
        event.setDesc(getNomenEventDesc());
        event.setOriginalRGDId(rgdId);
        event.setPreviousSymbol(gene.getSymbol());
        event.setPreviousName(gene.getName());
        dao.createNomenEvent(event);

        // change gene status to active
        rgd.setLastModifiedDate(new Date());
        rgd.setObjectStatus("ACTIVE");
        dao.updateRgdId(rgd);
        return true;
    }


    public void setFile(String file) {
        this.file = file;
    }

    /**
     * get a file with a list of RGD IDS to be reactivated;
     * sample line from the file:
     * <pre>
     * [2012-06-28 16:29:36,578] - RGDID=5045246, EGID=292228
     * [2012-06-28 16:29:39,212] - RGDID=2321568, EGID=100364190
     * </pre>
     * Note: the file is produced automatically by the pipeline in the log folder
     * @return
     */
    public String getFile() {
        return file;
    }

    public void setRefKey(String refKey) {
        this.refKey = refKey;
    }

    public String getRefKey() {
        return refKey;
    }

    public void setNomenEventType(String nomenEventType) {
        this.nomenEventType = nomenEventType;
    }

    public String getNomenEventType() {
        return nomenEventType;
    }

    public void setNomenEventDesc(String nomenEventDesc) {
        this.nomenEventDesc = nomenEventDesc;
    }

    public String getNomenEventDesc() {
        return nomenEventDesc;
    }
}
