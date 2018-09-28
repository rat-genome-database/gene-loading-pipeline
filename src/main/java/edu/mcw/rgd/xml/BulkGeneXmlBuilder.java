package edu.mcw.rgd.xml;

import edu.mcw.rgd.dataload.BulkGene;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.XdbManager;
import nu.xom.*;
import nu.xom.Node;
import org.jaxen.*;
import org.jaxen.xom.*;

import java.io.StringReader;
import java.sql.Timestamp;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: Mar 26, 2010
 * Time: 8:07:11 AM
 * it restores BulkGene object from the data serialized as XML
 */
public class BulkGeneXmlBuilder {

    public BulkGene run(String xml) throws Exception {

        BulkGene bg = new BulkGene();
        Gene gene = bg.getGene();

        // construct xom parser
        StringReader reader = new StringReader(xml);
        Builder builder = new Builder();
        Document doc = builder.build(reader);
        reader.close();
        Element root = doc.getRootElement();

        // recno
        bg.setRecNo(Integer.parseInt(root.getAttributeValue("recNo")));

        bg.setEgId(xpEGID.stringValueOf(root));
        bg.setEgType(xpEgType.stringValueOf(root));
        bg.geneOfficialSymbol = xpOfficialSymbol.stringValueOf(root);
        bg.geneOfficialName = xpOfficialName.stringValueOf(root);

        gene.setSpeciesTypeKey(SpeciesType.parse(xpSpecies.stringValueOf(root)));
        gene.setRgdId(xpRGDID.numberValueOf(root).intValue());
        gene.setSymbol(xpSymbol.stringValueOf(root));
        gene.setName(xpName.stringValueOf(root));
        gene.setDescription(xpDesc.stringValueOf(root));

        List<Node> aliases = xpAliases.selectNodes(root);
        for( Node node: aliases ) {
            Alias alias = new Alias();
            alias.setValue(node.getValue());
            alias.setTypeName("old_gene_symbol");
            alias.setNotes("From entrezgene " + (new Timestamp(System.currentTimeMillis())).toString());
            bg.aliases.addAlias(alias);
        }

        List<Element> maps = xpMaps.selectNodes(root);
        for( Element el: maps ) {

            MapData map = new MapData();

            String val = el.getAttributeValue("chromosome");
            map.setChromosome(val);
            if( val!=null && val.trim().length()> 0 )
                bg.setChromosome(val);

            val = el.getAttributeValue("fishBand");
            map.setFishBand(val);
            if( val!=null && val.trim().length()> 0 )
                bg.setFishband(map.getFishBand());

            val = el.getAttributeValue("startPos");
            if( val!=null && val.trim().length()>0 )
                map.setStartPos(Integer.parseInt(val));

            val = el.getAttributeValue("stopPos");
            if( val!=null && val.trim().length()>0 )
                map.setStopPos(Integer.parseInt(val));

            map.setStrand(el.getAttributeValue("strand"));
            map.setNotes("Crrated on " + (new Timestamp(System.currentTimeMillis())).toString());

            val = el.getAttributeValue("key");
            if( val!=null && val.trim().length()>0 )
                map.setKey(Integer.parseInt(val));

            val = el.getAttributeValue("srcPipeline");
            map.setSrcPipeline(val);

            bg.genePositions.addMapData(map);
        }

        List<Element> xdbs = xpXdbs.selectNodes(root);
        for( Element el: xdbs ) {

            List<XdbId> xdbids = bg.getXdbIds();
            XdbId xdb = new XdbId();
            xdb.setXdbKey(Integer.parseInt(el.getAttributeValue("key")));
            xdb.setAccId(el.getAttributeValue("id"));
            xdb.setSrcPipeline(XdbManager.EG_PIPELINE);
            xdb.setRgdId(gene.getRgdId());
            xdbids.add(xdb);
        }

        return bg;
    }

    static XPath xpEGID, xpSpecies, xpRGDID, xpEgType, xpSymbol, xpName, xpDesc;
    static XPath xpOfficialSymbol, xpOfficialName;
    static XPath xpAliases, xpMaps, xpXdbs;
    static {
        try {
            xpEGID = new XOMXPath("/EntrezGene/EntrezgeneID");
            xpSpecies = new XOMXPath("/EntrezGene/Species");
            xpRGDID = new XOMXPath("/EntrezGene/RGDID");
            xpEgType = new XOMXPath("/EntrezGene/EntrezgeneType");
            xpSymbol = new XOMXPath("/EntrezGene/GeneSymbol");
            xpName = new XOMXPath("/EntrezGene/GeneName");
            xpDesc = new XOMXPath("/EntrezGene/GeneDescription");
            xpOfficialSymbol = new XOMXPath("/EntrezGene/OfficialSymbol");
            xpOfficialName = new XOMXPath("/EntrezGene/OfficialName");

            xpAliases = new XOMXPath("/EntrezGene/Aliases/AliasValue");
            xpMaps = new XOMXPath("/EntrezGene/Maps/Map");
            xpXdbs = new XOMXPath("/EntrezGene/XDBSet/XDB");

        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
}
