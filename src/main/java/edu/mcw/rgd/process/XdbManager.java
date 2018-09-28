package edu.mcw.rgd.process;

import edu.mcw.rgd.datamodel.XdbId;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: Mar 3, 2010
 * Time: 10:18:44 AM
 * handles common logic when working with external ids;
 * f.e. EntrezGene ID and RGD ID
 */
public class XdbManager {

    public static final String EG_PIPELINE = "ENTREZGENE"; // name of entrezgene pipeline as found in RGD_ACC_XDB

    // remove GeneBank nucleotides (accId starts with AC_,NC_,NT_,NW_,NZ_);
    // never returns null -- just empty list
    public List<XdbId> removeGeneBankNucleotides(List<XdbId> ids) {

        if( ids==null )
            return new ArrayList<XdbId>();
        
        Iterator it = ids.listIterator();
        while( it.hasNext() ) {
            XdbId xdb = (XdbId) it.next();
            if( xdb.getXdbKey()!=XdbId.XDB_KEY_GENEBANKNU )
                continue;
            String accId = xdb.getAccId();
            if( accId==null )
                continue;
            if( accId.startsWith("AC_") || accId.startsWith("NC_") || accId.startsWith("NT_")
                    || accId.startsWith("NW_") || accId.startsWith("NZ_"))
            {
                it.remove();
            }
        }
        return ids;
    }
}
