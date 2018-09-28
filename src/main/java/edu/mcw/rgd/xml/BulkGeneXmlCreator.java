package edu.mcw.rgd.xml;

import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.dataload.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import nu.xom.*;
import org.apache.commons.logging.*;

import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 *
 * write BulkGene object in XML format
 *
 * uses singleton model to efficiently cache XDB names (to avoid numerous calls to db to resolve name from id)
 */
public class BulkGeneXmlCreator {

    private Document doc; // XOM document

    protected final Log logger = LogFactory.getLog("process");

    private XdbIdDAO xdbDAO = new XdbIdDAO();

    // internal map of XDB_ID to XDB_NAME
    private Map<Integer, String> _mapXdbIdToXdbName;

    // singleton access
    static private BulkGeneXmlCreator _instance = new BulkGeneXmlCreator();
    private BulkGeneXmlCreator() {
        _mapXdbIdToXdbName = new HashMap<Integer, String>();
    }
    static public BulkGeneXmlCreator getInstance() {
        return _instance;
    }

    // convert a single bulk gene to xml string
    public String createXmlString(BulkGene bg){
        Element bgEle = createBulkGeneElement(bg);
        return bgEle.toXML();
    }
    
    // convert bulk gene list to xml string
    public String createXmlString(ArrayList<BulkGene> bgs) {

        createDocument();
        createDOMTree(bgs);

        String xml = doc.toXML();
        logger.debug("Created xml String from bulk gene: " + xml);
        return xml;
    }
    
    private void createDocument() {
            
        Element root = new Element("EntrezGeneSet");
        root.addAttribute(new Attribute("ver", "1.0.3"));
        doc = new Document(root);
    }

    // create bulk gene dom tree for a list of bulk genes
    private void createDOMTree(ArrayList<BulkGene> bgs){
        //create the root element <EntrezGeneSet>
        Element rootEle = doc.getRootElement();

        //For each BulkGene object  create <Entrezgene> element and attach it to root
        for(BulkGene bg: bgs) {
            Element bgEle = createBulkGeneElement(bg);
            rootEle.appendChild(bgEle);
        }
    }   
    
    private Element createBulkGeneElement(BulkGene bg) {

        Gene gene = bg.getGene();
        Element bgEle = new Element("EntrezGene");
        bgEle.addAttribute(new Attribute("recNo", Integer.toString(bg.getRecNo())));

        //create EntrezgeneID element
        Element el = new Element("EntrezgeneID");
        el.appendChild(new Text(bg.getEgId()));
        bgEle.appendChild(el);

        //create Species element
        el = new Element("Species");
        el.appendChild(new Text(SpeciesType.getTaxonomicName(gene.getSpeciesTypeKey())));
        bgEle.appendChild(el);

        //create RGDID element
        el = new Element("RGDID");
        el.appendChild(new Text(Integer.toString(gene.getRgdId())));
        bgEle.appendChild(el);

        //create optional MatchingRGDID element
        if( bg.rgdGene!=null && bg.rgdGene.getRgdId()!=0 ) {
            el = new Element("MatchingRGDID");
            el.appendChild(new Text(Integer.toString(bg.rgdGene.getRgdId())));
            bgEle.appendChild(el);
        }

        //create EntrezGeneType element
        el = new Element("EntrezgeneType");
        el.appendChild(new Text(bg.getEgType()));
        bgEle.appendChild(el);

        //create RefSeq Status element
        if( gene.getRefSeqStatus()!=null ) {
            el = new Element("RefSeqStatus");
            el.appendChild(new Text(gene.getRefSeqStatus()));
            bgEle.appendChild(el);
        }

        //create NCBI Annot Status element
        if( gene.getNcbiAnnotStatus()!=null ) {
            el = new Element("NcbiAnnotStatus");
            el.appendChild(new Text(gene.getNcbiAnnotStatus()));
            bgEle.appendChild(el);
        }

        //create GeneSymbol element
        el = new Element("GeneSymbol");
        if( bg.rgdGene!=null && !Utils.stringsAreEqual(gene.getSymbol(), bg.rgdGene.getSymbol()) ) {
            el.addAttribute(new Attribute("Incoming", Utils.NVL(gene.getSymbol(),"")));
            el.addAttribute(new Attribute("InRgd", Utils.NVL(bg.rgdGene.getSymbol(),"")));
        } else {
            el.appendChild(new Text(gene.getSymbol()));
        }
        bgEle.appendChild(el);

        //create GeneName element
        el = new Element("GeneName");
        if( bg.rgdGene!=null && !Utils.stringsAreEqual(gene.getName(), bg.rgdGene.getName()) ) {
            el.addAttribute(new Attribute("Incoming", Utils.NVL(gene.getName(),"")));
            el.addAttribute(new Attribute("InRgd", Utils.NVL(bg.rgdGene.getName(),"")));
        } else {
            el.appendChild(new Text(gene.getName()));
        }
        bgEle.appendChild(el);

        //create GeneDesc element
        el = new Element("GeneDesc");
        if( bg.rgdGene!=null && !Utils.stringsAreEqual(gene.getDescription(), bg.rgdGene.getDescription()) ) {
            el.addAttribute(new Attribute("Incoming", Utils.NVL(gene.getDescription(),"")));

            // to avoid "nu.xom.IllegalCharacterDataException: 0x16 is not allowed in XML content"
            String desc = Utils.NVL(bg.rgdGene.getDescription(),"");
            if( desc.indexOf('\u0016')>=0 )
                desc = desc.replace('\u0016', ' ');

            el.addAttribute(new Attribute("InRgd", desc));
        } else {
            el.appendChild(new Text(gene.getDescription()));
        }
        bgEle.appendChild(el);

        // aliases
        el = new Element("Aliases");
        for( Alias alias: bg.aliases.getIncoming() ) {
            Element aEle = new Element("AliasValue");
            aEle.appendChild(new Text(alias.getValue()));
            aEle.addAttribute(new Attribute("type", alias.getTypeName()));
            el.appendChild(aEle);
        }
        bgEle.appendChild(el);

        // create map element and attach it to the bulkGeneElement
        if( bg.genePositions.isMapDataChanged() ) {
            el = new Element("Maps");
            bgEle.appendChild(el);

            createMaps("Incoming", bg.genePositions.getMapData(), el, false);
            createMaps("MapsInRgd", bg.genePositions.getMapDataInRgd(), el, true);
            createMaps("MapPosInserted", bg.genePositions.getMdForInsert(), el, false);
            createMaps("MapPosDeleted", bg.genePositions.getMdForDelete(), el, true);
        }
        else {
            createMaps("Maps", bg.genePositions.getMapData(), bgEle, false);
        }

        // create XDB element and attach it to the bulkGeneElement
        sortXdbIds(bg.getXdbIds()); // sort xdb-ids by XDB_KEY first, then by ACC_ID
        
        el = new Element("XDBSet");
        for(XdbId xdb: bg.getXdbIds()) {
            Element xdbEle = createXdbElement(xdb);
            el.appendChild(xdbEle);
        }
        bgEle.appendChild(el);

        el = new Element("Transcripts");
        for( TranscriptInfo tr: bg.transcripts ) {
            Element tEle = createTranscriptElement(tr);
            el.appendChild(tEle);
        }
        bgEle.appendChild(el);

        return bgEle;
    }

    private void createMaps(String elName, List<MapData> mds, Element parentElement, boolean exportMapKey) {
        if( mds==null || mds.isEmpty() )
            return;

        Element el = new Element(elName);
        for(MapData mapData: mds) {
            Element mapEle = createMapElement(mapData, exportMapKey);
            el.appendChild(mapEle);
        }
        parentElement.appendChild(el);
    }

    private Element createXdbElement(XdbId xdb) {

        Element el = new Element("XDB");

        //create xdb name element and text node and attach it to bulkGeneElement
        Attribute attr = new Attribute("key", Integer.toString(xdb.getXdbKey()));
        el.addAttribute(attr);

        //create xdb name element and text node and attach it to bulkGeneElement
        String xdbName = getXdbNameForXdbKey(xdb.getXdbKey());
        attr = new Attribute("name", xdbName);
        el.addAttribute(attr);

        //create xdb id element
        if( xdb.getAccId()!=null ) {
            attr = new Attribute("id", xdb.getAccId());
            el.addAttribute(attr);
        }
        return el;
    }
    
    private Element createMapElement(MapData map, boolean exportMapKey) {

        Element el = new Element("Map");

        //create map name element and text node and attach it to bulkGeneElement
        el.addAttribute(new Attribute("mapKey", Integer.toString(map.getMapKey())));
        if( map.getChromosome()!=null )
            el.addAttribute(new Attribute("chromosome", map.getChromosome()));
        if( map.getFishBand()!=null )
            el.addAttribute(new Attribute("fishBand", map.getFishBand()));
        if( map.getBandType()!=null )
            el.addAttribute(new Attribute("bandType", map.getBandType()));
        if( map.getStartPos()!=null )
            el.addAttribute(new Attribute("startPos", Integer.toString(map.getStartPos())));
        if( map.getStopPos()!=null )
            el.addAttribute(new Attribute("stopPos", Integer.toString(map.getStopPos())));
        if( map.getStrand()!=null )
            el.addAttribute(new Attribute("strand", map.getStrand()));
        if( map.getSrcPipeline()!=null )
            el.addAttribute(new Attribute("srcPipeline", map.getSrcPipeline()));
        if( exportMapKey )
            el.addAttribute(new Attribute("mapsDataKey", Integer.toString(map.getKey())));
        return el;
    }

    private Element createTranscriptElement(TranscriptInfo tr) {

        Element el = new Element("Tr");

        // transcript accession id
        Attribute attr = new Attribute("accId", tr.getAccId());
        el.addAttribute(attr);

        // transcript rgd id
        if( tr.getRgdId()>0 ) {
            attr = new Attribute("rgdId", Integer.toString(tr.getRgdId()));
            el.addAttribute(attr);
        }

        // protein accession id
        if( tr.getProteinAccId()!=null ) {
            attr = new Attribute("proteinAccId", tr.getProteinAccId());
            el.addAttribute(attr);
        }

        if( tr.getPeptideLabel()!=null && !tr.getPeptideLabel().isEmpty() ) {
            attr = new Attribute("peptideLabel", tr.getPeptideLabel());
            el.addAttribute(attr);
        }

        attr = new Attribute("coding", tr.isNonCoding()?"no":"yes");
        el.addAttribute(attr);
        if( tr.getRefSeqStatus()!=null ) {
            attr = new Attribute("refSeqStatus", tr.getRefSeqStatus());
            el.addAttribute(attr);
        }

        for( TranscriptLocus locus: tr.getLoci() ) {

            Element eloc = new Element("locus");
            el.appendChild(eloc);

            // transcript coords and map key
            attr = new Attribute("coords", dumpCoordsToString(locus.getTranscriptCoords()));
            eloc.addAttribute(attr);
            attr = new Attribute("mapKey", locus.getTranscriptCoords().getMapKey().toString());
            eloc.addAttribute(attr);
            attr = new Attribute("coords", dumpCoordsToString(locus.getTranscriptCoords()));
            eloc.addAttribute(attr);

            // dump utr5 coords
            if( locus.getUtr5()!=null ) {
                Element elc = new Element("utr5");
                attr = new Attribute("coords", dumpCoordsToString(locus.getUtr5()));
                elc.addAttribute(attr);
                eloc.appendChild(elc);
            }

            // dump utr3 coords
            if( locus.getUtr3()!=null ) {
                Element elc = new Element("utr3");
                attr = new Attribute("coords", dumpCoordsToString(locus.getUtr3()));
                elc.addAttribute(attr);
                eloc.appendChild(elc);
            }

            // dump exon coords
            for( MapData md: locus.getCoords() ) {
                Element elc = new Element("exon");
                attr = new Attribute("coords", dumpCoordsToString(md));
                elc.addAttribute(attr);
                eloc.appendChild(elc);
            }
        }
        return el;
    }

    String dumpCoordsToString(MapData md) {
        return md.getStartPos()+".."+md.getStopPos()+" ("+md.getStrand()+')';
    }

    // optimized for multiple calls
    String getXdbNameForXdbKey(int xdbKey) {

        // first check if the XDB_KEY is in the map
        String xdbName = _mapXdbIdToXdbName.get(xdbKey);
        if( xdbName!=null )
            return xdbName; // yes, it is -- we saved on expensive call to database

        // it is not -- get the value from database
        try {
            xdbName = xdbDAO.getXdbName(xdbKey);
            _mapXdbIdToXdbName.put(xdbKey, xdbName);
        }
        catch(Exception e) {
            e.printStackTrace();
            xdbName = ""; // exception ??
        }
        return xdbName;
    }

    // sort XdbIds by XDB_ID (or NAME) and then by acc_id
    void sortXdbIds(List<XdbId> xdbIds) {
        Collections.sort(xdbIds, new Comparator<XdbId>(){
            public int compare(XdbId x1, XdbId x2) {
                // first compare by xdb name
                int r = x1.getXdbKey() - x2.getXdbKey();
                if( r!=0 )
                    return r;
                r = Utils.stringsCompareTo(x1.getAccId(), x2.getAccId());
                if( r!=0 )
                    return r;
                return 0;
            }
        }
        );
    }

}
