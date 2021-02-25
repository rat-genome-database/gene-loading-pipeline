package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogFlagManager;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * @author mtutaj
 * @since 3/15/12
 * <p>
 * Handles everything related to positions of genes for one gene record. One gene record can have positions on multiple maps,
 * each identified by different MAP_KEY. Incoming positions are compared against existing positions in RGD and then
 * updated until they are in sync. Rule to avoid deleting of good data: only maps that appear in incoming data are synchronized,
 * to avoid deleting historical data. Also all database operations regarding gene positions are logged into gene_positions file
 */
public class GenePositions {
    private final Logger logger = Logger.getLogger("gene_positions");

    private List<MapData> mapData = new ArrayList<MapData>();
    private List<MapData> rgdMapData;

    EGDAO egDAO;

    public GenePositions(EGDAO dao) {
        egDAO = dao;
    }

    // list with results after running QC:
    // map positions matching perfectly
    // positions to be updated, to be inserted, and to be deleted
    List<MapData> mdMatching = new ArrayList<MapData>();
    List<MapData> mdForInsert = new ArrayList<MapData>();
    List<MapData> mdForDelete = new ArrayList<MapData>();

    /**
     * get incoming map data
     * @return
     */
    public List<MapData> getMapData() {
        return mapData;
    }

    /**
     * get map data in RGD
     * @return
     */
    public List<MapData> getMapDataInRgd() {
        return rgdMapData;
    }

    /**
     * check if incoming map data data is different than map data in RGD
     * @return true if incoming map data data is different than map data in RGD
     */
    public boolean isMapDataChanged() {
        return mdForInsert.size() + mdForDelete.size() > 0;
    }

    public List<MapData> getMdForInsert() {
        return mdForInsert;
    }

    public List<MapData> getMdForDelete() {
        return mdForDelete;
    }

    /**
     * add a new map data; a gene can have multiple positions on a map;
     * pseudoautosomal gens have positions on both chromosome X and Y;
     * many genes, due to missembly, have multiple positions on same chromosome
     * duplicate map positions won't be added
     * @param mdNew MapData to add
     * @return true if a position has been added; false if md was a duplicate
     */
    public boolean addMapData(MapData mdNew) {

        if( mdNew==null ) {
            System.out.println("null MapData");
            return false;
        }

        if( mdNew.getStartPos()!=null && mdNew.getStopPos()!=null && mdNew.getStartPos()==1 && mdNew.getStopPos()==1 ) {
            System.out.println("start pos==stop pos==1");
            return false;
        }

        for( MapData md: mapData ) {
            // look for duplicates: same map_key and chromosome and rgd_id
            if( md.equalsByGenomicCoords(mdNew) && md.getRgdId()==mdNew.getRgdId() ) {
                // duplicate found
                return false;
            }
        }

        // no duplicates -- add new position
        mapData.add(mdNew);
        return true;
    }

    public boolean hasGenomicPositionForAssembly(int mapKey) {
        for( MapData md: mapData ) {
            if( md.getMapKey()==mapKey
              && !Utils.isStringEmpty(md.getChromosome())
              && Utils.intsCompareTo(md.getStartPos(),0)>0
              && Utils.intsCompareTo(md.getStopPos(),0)>0 ) {
                return true;
            }
        }
        return false;
    }

    /**
     * update chromosome if not set
     * @param chr chromosome
     */
    public void updateChromosome(String chr) {

        for( MapData md: mapData ) {
            if( md.getChromosome()==null )
                md.setChromosome(chr);
        }
    }

    public void updateRgdId(int rgdId) {
        for( MapData md: mapData ) {
            md.setRgdId(rgdId);
        }
    }

    /**
     * load MapData for given RgdId;
     * remove from incoming data map keys not supported by the pipeline
     * @param rgdId rgd id
     * @param validMapKeys set of valid map keys
     */
    public void loadRgdMapData(int rgdId, Set<Integer> validMapKeys) throws Exception  {

        setRgdMapData(egDAO.getMapData(rgdId), validMapKeys);
    }

    /**
     * load MapData for given RgdId;
     * remove from incoming data map keys not supported by the pipeline
     * @param rgdMapData
     * @param validMapKeys set of valid map keys
     */
    public void setRgdMapData(List<MapData> rgdMapData, Set<Integer> validMapKeys) throws Exception  {

        this.rgdMapData = rgdMapData;
        if( validMapKeys!=null )
            restrictMapKeys(rgdMapData, validMapKeys);
    }

    /**
     * remove from incoming data map keys not supported by the pipeline
     * @param mdList - list of MapData object
     * @param validMapKeys set of valid map keys
     */
    private void restrictMapKeys(List<MapData> mdList, Set<Integer> validMapKeys) {

        // get a copy of rgd map data now, so data load thread will have the data ready
        // note: include only those map entries that are on official list of maps processed by pipeline
        Iterator<MapData> it = mdList.iterator();
        while( it.hasNext() ) {
            MapData md = it.next();
            boolean mapKeyMatches = validMapKeys.contains(md.getMapKey());
            if( !mapKeyMatches )
                // maps data is for assembly map not handled by the pipeline; remove it
                // so data loading module will not consider this map entry
                it.remove();
        }
    }

    /**
     * qc map data between rgd and incoming data; only maps that are present in incoming data are qc-ed
     */
    public void qcMapData(BulkGene bg, Logger logger) throws Exception {

        // logger override
        if( logger==null )
            logger = this.logger;

        // make a copy of both data in RGD and incoming data
        List<MapData> mdIncomingCopy = new ArrayList<MapData>(mapData);
        List<MapData> mdRgdCopy = rgdMapData!=null ? new ArrayList<MapData>(rgdMapData) : new ArrayList<MapData>();

        // initialize result collections
        mdMatching.clear();
        mdForDelete.clear();
        mdForInsert.clear();

        // synchronize every map present in incoming data
        while( !mdIncomingCopy.isEmpty() ) {
            // get map positions for one map
            int mapKey = mdIncomingCopy.get(0).getMapKey();

            List<MapData> mdsIncoming = detachMapData(mapKey, mdIncomingCopy);
            List<MapData> mdsRgd = detachMapData(mapKey, mdRgdCopy);
            List<MapData> mdsForInsert = new ArrayList<MapData>();

            // process all positions for incoming data
            for( MapData md: mdsIncoming ) {
                // check if position matches exactly by coordinates a position in RGD
                MapData mdMatch = detachMatchingMapData(md, mdsRgd);
                if( mdMatch!=null )
                    mdMatching.add(mdMatch);
                else {
                    if( md.getRgdId()==0 )
                        md.setRgdId(bg.getCustomFlags().getRgdId());
                    mdsForInsert.add(md);
                }
            }

            // whatever is left, goes into insert or delete arrays
            if( !mdsForInsert.isEmpty() )
                mdForInsert.addAll(mdsForInsert);
            if( !mdsRgd.isEmpty() )
                mdForDelete.addAll(mdsRgd);
        }
    }

    /**
     * synchronize map data between rgd and incoming data; only maps that are present in incoming data are synced
     */
    public void syncMapData(BulkGene bg, Logger logger, PipelineLogFlagManager dbFlagManager, String counterPrefix,
                            boolean keepOnePos, CounterPool counters) throws Exception {

        // logger override
        if( logger==null )
            logger = this.logger;

        // handle matching positions
        if( !mdMatching.isEmpty() ) {
            counters.add(counterPrefix+"_MAPPOS_MATCHED", mdMatching.size());
        }

        // handle position updates
        if( !mdForInsert.isEmpty() ) {

            egDAO.insertMapData(mdForInsert);

            for( MapData md: mdForInsert ) {
                logger.info("EGID="+bg.getEgId()+" MAPS_DATA INSERT >"+dumpMapPosition(md));
            }

            dbFlagManager.setFlag(counterPrefix+"_MAPPOS_INSERTED", bg.getRecNo());
            counters.add(counterPrefix+"_MAPPOS_INSERTED", mdForInsert.size());
        }

        // handle position deletes
        if( !mdForDelete.isEmpty() ) {

            if( keepOnePos ) {
                int deletionsSuppressed = suppressPositionsForDelete();
                if( deletionsSuppressed>0 ) {
                    dbFlagManager.setFlag(counterPrefix + "_MAPPOS_DELETE_SUPPRESSED", bg.getRecNo());
                    counters.add(counterPrefix + "_MAPPOS_DELETE_SUPPRESSED", deletionsSuppressed);
                }
            }

            if( egDAO.deleteMapData(mdForDelete, counterPrefix.equals("GENE")) != 0 ) {

                for (MapData md : mdForDelete) {
                    logger.info("EGID=" + bg.getEgId() + " MAPS_DATA DELETE >" + dumpMapPosition(md));
                }

                dbFlagManager.setFlag(counterPrefix + "_MAPPOS_DELETED", bg.getRecNo());
                counters.add(counterPrefix + "_MAPPOS_DELETED", mdForDelete.size());
            } else {
                dbFlagManager.setFlag(counterPrefix + "_MAPPOS_DELETE_SUPPRESSED", bg.getRecNo());
                counters.add(counterPrefix + "_MAPPOS_DELETE_SUPPRESSED", mdForDelete.size());

            }
        }
    }

    private int suppressPositionsForDelete() {
        int suppressed = 0;
        Iterator<MapData> it = mdForDelete.iterator();
        while( it.hasNext() ) {
            MapData md = it.next();

            // check if there is a similar position for same rgd_id, map_key
            // among matches and inserts
            if( !findSimilarPosition(md, mdMatching) &&
                    !findSimilarPosition(md, mdForInsert) ) {
                logger.info(" MAPS_DATA DELETE_SUPPRESSED >"+dumpMapPosition(md));
                suppressed++;
                it.remove();
            }
        }
        return suppressed;
    }

    private boolean findSimilarPosition(MapData md, List<MapData> mds) {

        for( MapData md2: mds ) {
            if( md.getRgdId()==md2.getRgdId() && md.getMapKey().equals(md2.getMapKey()) ) {
                return true;
            }
        }
        return false;
    }

    private String dumpMapPosition(MapData md) {

        return "MAPS_DATA_KEY:"+md.getKey()+
              "|CHROMOSOME:"+md.getChromosome()+
              "|MAP_KEY:"+md.getMapKey()+
              "|RGD_ID:"+md.getRgdId()+
              "|START_POS:"+md.getStartPos()+
              "|STOP_POS:"+md.getStopPos()+
              "|STRAND:"+md.getStrand()+
              "|SRC_PIPELINE:"+md.getSrcPipeline()+
              "|FISH_BAND:"+md.getFishBand()+
              "|ABS_POSITION:"+md.getAbsPosition()+
              "|NOTES:"+md.getNotes()+
            "";
    }

    /**
     * detach from list of MapData object a subset matching given map key
     * @param mapKey map key
     * @param mds List of MapData objects
     * @return List<MapData>
     */
    private List<MapData> detachMapData(int mapKey, List<MapData> mds) {

        List<MapData> results = new ArrayList<MapData>(mds.size());

        Iterator<MapData> it = mds.iterator();
        while( it.hasNext() ) {
            MapData md = it.next();
            if( md.getMapKey()==mapKey ) {
                it.remove();
                results.add(md);
            }
        }
        return results;
    }

    public void deleteOverlappingPositionsMarkedForDelete(BulkGene bg, PipelineLogFlagManager dbFlagManager, CounterPool counters) throws Exception {

        List<MapData> mdsOverlapping = new ArrayList<MapData>();

        for( MapData mdForDelete: this.getMdForDelete() ) {

            // find a matching overlapping position in-rgd
            for( MapData mdInRgd: this.getMapDataInRgd() ) {
                // skip self-matches
                if( mdForDelete.getKey()==mdInRgd.getKey() )
                    continue;
                if( positionsOverlap(mdForDelete, mdInRgd) ) {
                    // for-delete position overlaps with another in-rgd position
                    // this is most likely a legitimate update of transcript position
                    // so delete the old position
                    mdsOverlapping.add(mdForDelete);
                    break;
                }
            }
        }

        if( !mdsOverlapping.isEmpty() ) {
            if( egDAO.deleteMapData(mdsOverlapping, false)!=0 ) {
                for( MapData md: mdsOverlapping ) {
                    logger.info("OVERLAPPING MAPS_DATA DELETE >"+dumpMapPosition(md));
                }

                dbFlagManager.setFlag("TRANSCRIPT_OVERLAPPING_MAPPOS_DELETED", bg.getRecNo());
                counters.add("TRANSCRIPT_OVERLAPPING_MAPPOS_DELETED", mdsOverlapping.size());
            } else {
                dbFlagManager.setFlag("TRANSCRIPT_OVERLAPPING_MAPPOS_DELETE_SUPPRESSED", bg.getRecNo());
                counters.add("TRANSCRIPT_OVERLAPPING_MAPPOS_DELETE_SUPPRESSED", mdsOverlapping.size());
            }

        }
    }

    public static boolean positionsOverlap(MapData md1, MapData md2) {

        if(!md1.getMapKey().equals(md2.getMapKey()))
            return false;
        if( !md1.getChromosome().equals(md2.getChromosome()) )
             return false;
        if( md1.getStopPos() < md2.getStartPos() )
            return false;
        if( md2.getStopPos() < md1.getStartPos() )
            return false;
        return true;
    }

    MapData getPositionMatchingLocus(MapData locus) {
        for( MapData md: mapData ) {
            if( positionsOverlap(md, locus) ) {
                return md;
            }
        }
        return null;
    }

    /**
     * detach from list of MapData object a MapData object matching by coordinates and rgdId
     * @param md MapData
     * @param mds List of MapData objects
     * @return List<MapData>
     */
    private MapData detachMatchingMapData(MapData md, List<MapData> mds) {

        Iterator<MapData> it = mds.iterator();
        while( it.hasNext() ) {
            MapData md2 = it.next();
            if( md.equalsByGenomicCoords(md2) ) {

                // match by position; subject to detach unless md and md2 has non-zero rgdIds that are different
                if( md.getRgdId()<=0 || md2.getRgdId()<=0 || md.getRgdId()==md2.getRgdId() ) {
                    it.remove();
                    return md2;
                }
            }
        }
        return null;
    }
}
