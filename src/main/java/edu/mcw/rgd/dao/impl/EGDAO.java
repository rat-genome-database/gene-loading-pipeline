package edu.mcw.rgd.dao.impl;

import edu.mcw.rgd.dao.spring.IntListQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.Map;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.*;
import java.util.*;

/**
 * @author mtutaj
 * @since Mar 8, 2010
 * DAO code specific to the pipeline (singleton)
 */
public class EGDAO {

    private AliasDAO aliasDAO = new AliasDAO();
    private AssociationDAO assocDAO = new AssociationDAO();
    private GeneDAO geneDAO = assocDAO.getGeneDAO();
    private MapDAO mapDAO = new MapDAO();
    private NomenclatureDAO nomenclatureDAO = new NomenclatureDAO();
    private RGDManagementDAO rgdDAO = new RGDManagementDAO();
    private TranscriptDAO transcriptDAO = new TranscriptDAO();
    private XdbIdDAO xdbidDAO = new XdbIdDAO();

    protected final Logger logXdbIds = Logger.getLogger("xdb_ids");
    protected final Logger logAliases = Logger.getLogger("aliases");
    protected final Logger logAssoc = Logger.getLogger("assoc");

    static EGDAO _instance;

    static public EGDAO getInstance() {
        if( _instance==null ) {
            _instance = new EGDAO();
        }
        return _instance;
    }

    private EGDAO() {
    }

    public List<Integer> getEgIdsForGenesWithoutStrand(int speciesTypeKey, int mapKey) throws Exception {
        String sql = "SELECT DISTINCT TO_NUMBER(acc_id) "+
            "FROM rgd_acc_xdb x,rgd_ids r,maps_data md "+
            "WHERE x.xdb_key=3 AND x.rgd_id=r.rgd_id AND r.object_key=1 AND r.species_type_key=? AND r.object_status='ACTIVE' "+
            "AND md.rgd_id=r.rgd_id AND md.map_key=? AND md.strand IS NULL";
        return IntListQuery.execute(geneDAO, sql, speciesTypeKey, mapKey);
    }

    // get gene by symbol and species; return null if not found
    public Gene getGeneBySymbolAndSpecies(String geneSymbol, int speciesKey) throws Exception {

        List<Gene> genes = getGenesBySymbolAndSpecies(geneSymbol, speciesKey);
        if( genes!=null && genes.size()>0 )
            return genes.get(0);
        return null;
    }

    // get gene by symbol and species; return null if not found
    public List<Gene> getGenesBySymbolAndSpecies(String geneSymbol, int speciesKey) throws Exception {
        return geneDAO.getAllGenesBySymbol(geneSymbol, speciesKey);
    }

    // get all RGD_ACC_XDB ids records from database; the fields will contain:
    // ACC_ID, XDB_KEY, RGD_ID, OBJECT_STATUS, SPECIES_TYPE_KEY
    public ResultSet getAllXdbIds() throws Exception {

        String query = "select a.ACC_ID,a.XDB_KEY,a.RGD_ID,r.OBJECT_STATUS,r.SPECIES_TYPE_KEY from RGD_ACC_XDB a,RGD_IDS r "+
                "where a.RGD_ID=r.RGD_ID and a.ACC_ID like '__!_%' escape '!'";
        
        Connection conn = mapDAO.getConnection();
        PreparedStatement ps = conn.prepareStatement(query);
        return ps.executeQuery();
    }

    public int getPrimaryMapKey(int speciesTypeKey) throws Exception {
        if( speciesTypeKey==8 ) {
            return 0; // no primary ref assembly for zebrafish
        }

        Map map = mapDAO.getPrimaryRefAssembly(speciesTypeKey);
        if( map!=null )
            return map.getKey();
        else
            return 0;
    }

    /**
     * get positional information for object identified by rgd id
     * @param rgdId rgd id of object
     * @return list of mapData objects
     * @throws Exception
     */
    public List<MapData> getMapData(int rgdId) throws Exception {
       return mapDAO.getMapData(rgdId);
    }

    /**
     * get positional information for object identified by rgd id and map_key
     * @param egId Ncbi GeneId
     * @param mapKey map key
     * @return count of MapData objects
     * @throws Exception
     */
    public int getMapDataCount(String egId, int mapKey) throws Exception {
        int count = 0;
        for( Integer rgdId: getRgdIdListByEGID(Integer.parseInt(egId)) ) {
            count += mapDAO.getMapData(rgdId, mapKey).size();
        }
        return count;
    }

    /**
     * insert a list of MapData objects
     * @param mdList list of MapData objects to be inserted
     * @return number of rows affected
     * @throws Exception if something unexpected happens in the framework
     */
    public int insertMapData(List<MapData> mdList) throws Exception{
        // always set src pipeline to 'NCBI'
        for( MapData md: mdList ) {
            md.setSrcPipeline("NCBI");
            // throw exception if map_key or rgd_id is not set
            if( md.getMapKey()==null || md.getMapKey()==0 || md.getRgdId()<=0 )
                throw new Exception("insert map data: no map key or no rgd id");
        }

        return mapDAO.insertMapData(mdList);

    }

    /**
     * delete map data given a list of map data key
     * @param mapDataList list of map data objects
     * @return number of rows deleted
     * @throws Exception if something wrong happens in spring framework
     */
    public int deleteMapData(List<MapData> mapDataList) throws Exception{
        return mapDAO.deleteMapData(mapDataList);
    }

    /**
     * get list of aliases for given RGD_ID; exclude 'array_id_' aliases
     * @param rgdId RGD_ID
     * @return list of aliases associated with given RGD_ID
     * @throws Exception if something wrong happens in spring framework
     */
    public List<Alias> getAliases(int rgdId) throws Exception {
        List<Alias> aliases = aliasDAO.getAliases(rgdId);
        Iterator<Alias> it = aliases.iterator();
        while( it.hasNext() ) {
            Alias alias = it.next();
            if( Utils.defaultString(alias.getTypeName()).startsWith("array_id") )
                it.remove();
        }
        return aliases;
    }

    /**
     * insert alias list into ALIASES table
     * @param aliases list of aliases to be inserted
     * @return count of rows affected
     * @throws Exception if something wrong happens in spring framework
     */
    public int insertAliases(List<Alias> aliases) throws Exception {
        if( aliases.isEmpty() ) {
            return 0;
        }

        for( Alias alias: aliases ) {
            logAliases.debug("INSERTING "+alias.dump("|"));
        }

        return aliasDAO.insertAliases(aliases);
    }

    /**
     * update last modified date for specified rgd id
     * @param rgdId rgd id
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void updateLastModifiedDate(int rgdId) throws Exception {
        rgdDAO.updateLastModifiedDate(rgdId);
    }

    /**
     * Update properties of rgd id object
     *
     * @param rgdId  rgd id
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void updateRgdId(RgdId rgdId) throws Exception {
        rgdDAO.updateRgdId(rgdId);
    }

    /**
     * create a new rgd_id object
     * @param objectKey object key
     * @param objectStatus object status: 'ACTIVE', 'RETIRED', 'WITHDRAWN'
     * @param notes additional notes about the application creating this rgd id
     * @param speciesTypeKey species type key: 1, 2, 3
     * @return newly created RgdId object; generated unique rgd_id value is returned
     * @throws Exception when unexpected error in spring framework occurs
     */
    public RgdId createRgdId(int objectKey, String objectStatus, String notes, int speciesTypeKey) throws Exception {
        return rgdDAO.createRgdId(objectKey, objectStatus, notes, speciesTypeKey);
    }

    /**
     * unlink feature from transcripts; the feature itself is not deleted
     * @param featureRgdId feature rgd id
     * @return number of rows affected
     * @throws Exception on error in framework
     */
    public int unlinkFeature(int featureRgdId) throws Exception {

        // delete the transcript feature itself
        String query = "delete from TRANSCRIPT_FEATURES where FEATURE_RGD_ID=?";
        return transcriptDAO.update(query, featureRgdId);
    }

    private List<Integer> getRgdIdListByEGID(int egId) throws Exception {

        List<Gene> genes = xdbidDAO.getGenesByXdbId(XdbId.XDB_KEY_ENTREZGENE, Integer.toString(egId));
        List<Integer> rgdIds = new ArrayList<>(genes.size());
        for( Gene gene: genes ) {
            // exclude splices and alleles from the results
            if( !gene.isVariant() )
                rgdIds.add(gene.getRgdId());
        }
        return rgdIds;
    }

    /// same as XdbIdDAO.getGenesByXdbId(), but all genes of type 'allele' or 'splice' are excluded
    public List<Gene> getGenesByEGID(String egId) throws Exception {

        List<Gene> genes = xdbidDAO.getGenesByXdbId(XdbId.XDB_KEY_ENTREZGENE, egId);
        Iterator<Gene> it = genes.listIterator();
        while( it.hasNext() ) {
            Gene gene = it.next();
            // remove all genes of type "splice" or "allele" from the list
            if( gene.isVariant() )
                it.remove();
        }
        return genes;
    }

    /**
     * get all genes with given external id
     * @param xdbKey external db key
     * @param accId external id to be looked for
     * @return list of Gene objects
     */
    public List<RgdId> getRGDIdsByXdbId(int xdbKey, String accId) throws Exception {

        return xdbidDAO.getRGDIdsByXdbId(xdbKey, accId);
    }

    public List<XdbId> getXdbIdsByRgdId(int xdbKey, int rgdId) throws Exception {

        return xdbidDAO.getXdbIdsByRgdId(xdbKey, rgdId);
    }

    /**
     * return external ids for any combination of parameters;
     * if given parameter is null or 0, it means, that any value of this parameter could be accepted
     *
     * @param filter acc_id,xdb_id,rgd_id and src_pipeline are checked
     * @return list of external ids
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<XdbId> getXdbIds(XdbId filter) throws Exception {

        return xdbidDAO.getXdbIds(filter);
    }

    /**
     * return external ids for any combination of parameters given in filter;
     * if given parameter is null or 0, it means, that any value of this parameter could be accepted
     *
     * @param filter any combination of acc_id,xdb_id,rgd_id and src_pipeline is honored
     * @param speciesType species type key
     * @param objectKey object key
     * @return list of external ids matching the filter; empty list is returned if no matching entries are found
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<XdbId> getXdbIds(XdbId filter, int speciesType, int objectKey) throws Exception {

        return xdbidDAO.getXdbIds(filter, speciesType, objectKey);
    }

    /**
     * delete a list external ids (RGD_ACC_XDB rows);
     * if ACC_XDB_KEY is provided, it is used to delete the row;
     * else ACC_ID, RGD_ID, XDB_KEY and SRC_PIPELINE are used to locate and delete every row
     *
     * @param xdbIds list of external ids to be deleted
     * @param objectType object type like 'TRANSCRIPT' or 'GENE' -- used in logging
     * @return nr of rows deleted
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int deleteXdbIds( Collection<XdbId> xdbIds, String objectType ) throws Exception {

        // sanity check
        if( xdbIds==null )
            return 0;

        for( XdbId xdbId: xdbIds ) {
            logXdbIds.info("DELETE "+objectType+"|"+xdbId.dump("|"));
        }
        return xdbidDAO.deleteXdbIds((List<XdbId>)xdbIds);
    }

    /**
     * update properties of list of XdbIds objects
     * @param xdbIds list of objects with data to be updated
     * @return count of rows affected
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int updateXdbIds(List<XdbId> xdbIds) throws Exception {

        // sanity check
        if( xdbIds==null )
            return 0;

        int affectedRows = 0;
        for( XdbId xdbId: xdbIds ) {
            logXdbIds.info("UPDATE|"+xdbId.dump("|"));
            affectedRows += xdbidDAO.updateByKey(xdbId);
        }
        return affectedRows;
    }

    /**
     * for a bunch of rows identified by acc_xdb_key, set MODIFICATION_DATE to SYSDATE
     * @return number of actually updated rows
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int updateModificationDate(Collection<XdbId> xdbIds) throws Exception {

        if( xdbIds==null || xdbIds.isEmpty() )
            return 0;

        List<Integer> accXdbKeys = new ArrayList<>(xdbIds.size());
        for( XdbId xdbId: xdbIds ) {
            accXdbKeys.add(xdbId.getKey());
        }
        return xdbidDAO.updateModificationDate(accXdbKeys);
    }

    /**
     * insert a bunch of XdbIds; duplicate entries are not inserted (with same RGD_ID,XDB_KEY,ACC_ID,SRC_PIPELINE)
     * @param xdbIds list of XdbIds objects to be inserted
     * @param objectType object type like 'TRANSCRIPT' or 'GENE' -- used in logging
     * @return number of actually inserted rows
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int insertXdbs(Collection<XdbId> xdbIds, String objectType) throws Exception {

        for( XdbId xdbId: xdbIds ) {
            logXdbIds.info("INSERT "+objectType+"|" + xdbId.dump("|"));
        }
        return xdbidDAO.insertXdbs((List<XdbId>)xdbIds);
    }

    /**
     * get active genes with given external id
     * @param xdbKey - external db key
     * @param accId - external id to be looked for
     * @return list of Gene objects
     */
    public List<Gene> getActiveGenesByXdbId(int xdbKey, String accId) throws Exception {

        return xdbidDAO.getActiveGenesByXdbId(xdbKey, accId);
    }

    /**
     * get RgdId object for given rgd id
     * @param rgdId rgd id of object
     * @return RgdId object
     * @throws Exception
     */
    public RgdId getRgdId(int rgdId) throws Exception {

        return rgdDAO.getRgdId2(rgdId);
    }

    /**
     * Returns a Gene based on an rgd id
     * @param rgdId rgd id
     * @return Gene object for given rgd id
     * @throws Exception thrown when there is no gene with such rgd id
     */
    public Gene getGene(int rgdId) throws Exception {

        return geneDAO.getGene(rgdId);
    }

    /**
     * Update gene in the datastore based on rgdID
     *
     * @param gene Gene object
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void updateGene(Gene gene) throws Exception{

        geneDAO.updateGene(gene);
    }

    /**
     * insert a new gene into GENES table; upon successful insert, 'gene' object will have 'key'
     * property set to the key assigned automatically
     * @param gene Gene object to be inserted
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void insertGene(Gene gene) throws Exception{

        geneDAO.insertGene(gene);
    }

    /**
     * insert new gene type into GENE_TYPES table
     * @param geneType new gene type
     * @param geneTypeDesc gene type description
     * @param geneDescPublic gene public description
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void createGeneType(String geneType, String geneTypeDesc, String geneDescPublic) throws Exception {

        geneDAO.createGeneType(geneType, geneTypeDesc, geneDescPublic);
    }

    // return true if given geneType is in GENE_TYPES table
    public boolean existsGeneType(String geneType) throws Exception {

        return geneDAO.existsGeneType(geneType);
    }

    public int getReplacedGeneByRgdId(int oldRgdId) throws Exception {
        return rgdDAO.getActiveRgdIdFromHistory(oldRgdId);
    }

    /**
     * given Entrez Gene Id, try to get corresponding Rgd id  (logic used by ortholog relation loading pipeline)
     * @param egId Entrez Gene Id
     * @return matching rgd id, <br>
     * or -1 unmatched ----- entrezgene id doesn't match any genes in RGD (including active and non active genes) <br>
     * or -2 multiple ---- entrezgene id matches multiple genes(excluding splices and alleles) or matches multiple genes that are all splice or allele genes <br>
     * or -3 withdrawn ---- entrezgene id matches non active gene which is not replaced by any active gene
     * @throws Exception
     */
    public int matchRgdIdByEgId(int egId) throws Exception {

        // the rgdIds queried by the entrezgene id; splice and allele are excluded
        List<Integer> rgdIds = getRgdIdListByEGID(egId);
        //logger.debug("rgd ids: "+ rgdIds);
        List<Integer> nonActiveRgdIds=new ArrayList<> ();

        if (rgdIds==null || rgdIds.isEmpty() )
            return -1; // unmatched

        String rRgdId=""; // the returned rgd id
        String rgdIdInfo=""; // the associated rgd ids if there is multiple match

        // the result has rgd id
        int activeC=0;     // count active genes
        for (int rgdId: rgdIds ) {
            // get status for each gene
            String status = getRgdId(rgdId).getObjectStatus();
            if (!status.equals("ACTIVE")) {
                nonActiveRgdIds.add(rgdId);
                continue;
            }
            rRgdId=rgdId+"";
            activeC++;
            rgdIdInfo=rgdIdInfo+","+rRgdId;
        }

        if ( activeC >1) {
            // the gene is associated with multiple active genes
            return -2; // multiple+":"+rgdIdInfo;
        } else if (activeC ==1) {
            // the gene matches only one active gene
            return Integer.parseInt(rRgdId);
        } else {
            // activeC=0, the gene is associated with only nonactive genes
            int replacedRgdId=0;
            //int replacedActiveC=0;
            String reRgdIdInfo="";

            // find the replaced genes
            List<Integer> replacedRgdIds= new ArrayList<> ();
            for (int nonActiveRgdId: nonActiveRgdIds ) {
                replacedRgdId=getReplacedGeneByRgdId(nonActiveRgdId);
                if (replacedRgdId >0) {
                    replacedRgdIds.add(replacedRgdId); // the replaced gene is active
                    //replacedActiveC++;
                    reRgdIdInfo=reRgdIdInfo+","+replacedRgdId;
                }
            }

            if (replacedRgdIds.size() >1) {
                    // multiple replaced genes
                    return -2; // multiple+" replaced:"+ reRgdIdInfo;
            } else if (replacedRgdIds.size() ==1) {
                // matches one replaced gene
                    return replacedRgdId;
            } else {
                // filteredReplacedRgdIds.size()=0, couldn't find replaced gene
                return -3; // withdrawn+":"+nonActiveRgdId;
            }
        }
    }

    public List<Association> getAssociationsByType(String assocType, int speciesTypeKey) throws Exception {
        return assocDAO.getAssociationsByType(assocType, speciesTypeKey);
    }

    public List<Association> getAssociationsByRgdId(int masterRgdId) throws Exception {
        return assocDAO.getAssociationsForMasterRgdId(masterRgdId);
    }

    public int insertAssociation(Association assoc) throws Exception {
        logAssoc.info("INSERT "+assoc.dump("|"));
        return assocDAO.insertAssociation(assoc);
    }

    public int updateAssociation(Association assoc) throws Exception {
        logAssoc.info("UPDATE "+assoc.dump("|"));
        return assocDAO.updateAssociation(assoc);
    }

    /**
     * delete an association given assoc key
     * @param assoc Association object to be deleted
     * @return count of rows affected: 1 on successful delete, 0 if invalid assoc key
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int deleteAssociation( Association assoc ) throws Exception {
        logAssoc.info("DELETE "+assoc.dump("|"));
        return assocDAO.deleteAssociationByKey(assoc.getAssocKey());
    }

    /**
     *  Creates a new nomenclature event in the datastore.
     * @param event NomenclatureEvent object
     * @throws Exception if something wrong happens in spring framework
     */
    public void createNomenEvent(NomenclatureEvent event) throws Exception {

        nomenclatureDAO.createNomenEvent(event);
    }

    public Collection<Integer> getEgIdsForAllActiveGenes(int speciesTypeKey) throws Exception {
        // get all EntrezGene ids for current species
        XdbId filter = new XdbId();
        filter.setXdbKey(XdbId.XDB_KEY_ENTREZGENE);
        EGDAO dao = EGDAO.getInstance();
        Set<Integer> egIds = new HashSet<>();
        for( XdbId xdbId: dao.getXdbIds(filter, speciesTypeKey, RgdId.OBJECT_KEY_GENES) ) {

            // skip GeneIds linked to inactive RGD objects
            RgdId rgdId = dao.getRgdId(xdbId.getRgdId());
            if( !rgdId.getObjectStatus().equals("ACTIVE") ) {
                continue;
            }

            egIds.add(Integer.parseInt(xdbId.getAccId()));
        }
        return egIds;
    }

    public Collection<Integer> getEgIdsForActiveGenesFromFile(int speciesTypeKey, String fileName) throws Exception {

        EGDAO dao = EGDAO.getInstance();
        Set<Integer> egIds = new HashSet<>();

        // load eg-ids from file
        String fileContent = Utils.readFileAsString(fileName);
        for( String line: fileContent.split("[\\r\\n]") ) {
            String egId = line.trim();
            int geneId = Integer.parseInt(egId);
            for( Gene gene: dao.getGenesByEGID(egId) ) {
                // skip different species
                if( gene.getSpeciesTypeKey()!=speciesTypeKey ) {
                    continue;
                }

                // skip GeneIds linked to inactive RGD objects
                RgdId rgdId = dao.getRgdId(gene.getRgdId());
                if( !rgdId.getObjectStatus().equals("ACTIVE") ) {
                    continue;
                }

                egIds.add(geneId);
                break;
            }
        }
        return egIds;
    }

    public Collection<Integer> getEdIdsForMicroRnaGenes(int speciesTypeKey) throws Exception {

        String sql = "SELECT DISTINCT TO_NUMBER(acc_id) "+
            "FROM rgd_acc_xdb x,rgd_ids r,genes g "+
            "WHERE x.xdb_key=3 AND x.rgd_id=r.rgd_id AND r.object_key=1 AND r.species_type_key=? AND r.object_status='ACTIVE' "+
            "AND g.rgd_id=r.rgd_id AND g.gene_symbol_lc LIKE 'mir%'";
        return IntListQuery.execute(geneDAO, sql, speciesTypeKey);
    }

    public Collection<Integer> getEdIdsForTRnaGenes(int speciesTypeKey) throws Exception {

        String sql = "SELECT DISTINCT TO_NUMBER(acc_id) "+
            "FROM rgd_acc_xdb x,rgd_ids r,genes g "+
            "WHERE x.xdb_key=3 AND x.rgd_id=r.rgd_id AND r.object_key=1 AND r.species_type_key=? AND r.object_status='ACTIVE' "+
            "AND g.rgd_id=r.rgd_id AND g.gene_type_lc='trna'";
        return IntListQuery.execute(geneDAO, sql, speciesTypeKey);
    }

    public Collection<Integer> getEdIdsForMergedGenes() throws Exception {

        ArrayList<Integer> egIds = new ArrayList<>(2000);
        BufferedReader reader = new BufferedReader(new FileReader("/tmp/rat_merges_7-24-15.txt"));
        String line = reader.readLine(); // skip header
        while( (line=reader.readLine())!=null ) {
            String[] cols = line.split("[\\t]",-1);
            Integer egId = Integer.parseInt(cols[2]);
            egIds.add(egId);
        }
        reader.close();
        return egIds;
    }

    /**
     * detach a transcript from gene, by removing a row from TRANSCRIPTS table
     * @param tr Transcript object
     * @return number of rows affected
     * @throws Exception on error in framework
     */
    public int detachTranscriptFromGene(Transcript tr) throws Exception {
        return transcriptDAO.detachTranscriptFromGene(tr.getRgdId(), tr.getGeneRgdId());
    }

    /**
     * update transcript object; TRANSCRIPT_RGD_ID will be used to identify the transcript
     * @param tr transcript object
     * @return number of rows affected
     * @throws Exception on error in spring framework
     */
    public int updateTranscript(Transcript tr) throws Exception {
        return transcriptDAO.updateTranscript(tr);
    }
    /**
     * create a new transcript object and assign genomic coordinates
     * @param tr transcript object
     * @param speciesTypeKey species type key
     * @return number of rows affected
     * @throws Exception on error in spring framework
     */
    public int createTranscript(Transcript tr, int speciesTypeKey) throws Exception {
        return transcriptDAO.createTranscript(tr, speciesTypeKey);
    }

    /**
     * bind existing feature object to a transcript
     * @param transcriptRgdId transcript rgd id
     * @param featureRgdId rgd id of transcript feature
     * @throws Exception on error in framework
     */
    public void bindFeatureToTranscript(int transcriptRgdId, int featureRgdId) throws Exception {
        transcriptDAO.bindFeatureToTranscript(transcriptRgdId, featureRgdId);
    }

    /**
     * create a new transcript feature object
     * @param tf transcript feature
     * @param speciesTypeKey species type key
     * @return number of rows affected
     * @throws Exception on error in framework
     */
    public int createFeature(TranscriptFeature tf, int speciesTypeKey) throws Exception {
        return transcriptDAO.createFeature(tf, speciesTypeKey);
    }
}
