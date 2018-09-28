package edu.mcw.rgd.process;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 7/27/15
 * Time: 10:10 AM
 * <p>
 * reads a list of genes to be merged from an external file, and then runs a merge
 * input file format:
 * merge-from-rgd-id, merge-to-rgd-id, keep-this-ncbi-gene-id
 */
public class BulkGeneMerge {

    /**
     *
     * @param args
     * @throws Exception
     * <p>
     * input file format:
    [2017-07-13 14:26:57,736] - GeneTrackStatus=SECONDARY|OldGeneId=101967056|GeneRGDId=12732137|GeneSymbol=LOC101967056|CurrentGeneId=101967252|Species=Squirrel     */
    public static void main(String[] args) throws Exception {

        EGDAO dao = EGDAO.getInstance();

        Map<Integer,Integer> mergeMap = new HashMap<>();
        int speciesTypeKey = 0;

        BufferedWriter out = new BufferedWriter(new FileWriter("/tmp/notmerged.txt"));
        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line;
        while( (line=reader.readLine())!=null ) {
            String[] cols = line.split("[\\|]", -1);

            String fromGeneId = null;
            if( cols[1].startsWith("OldGeneId=") ) {
                fromGeneId = cols[1].substring(10);
            }
            int fromRgdIdIncoming = 0;
            if( cols[2].startsWith("GeneRGDId=") ) {
                fromRgdIdIncoming = Integer.parseInt(cols[2].substring(10));
            }
            String fromGeneSymbol = null;
            if( cols[3].startsWith("GeneSymbol=") ) {
                fromGeneSymbol = cols[3].substring(11);
            }
            String toGeneId = null;
            if( cols[4].startsWith("CurrentGeneId=") ) {
                toGeneId = cols[4].substring(14);
            }
            speciesTypeKey = SpeciesType.ALL;
            if( cols[5].startsWith("Species=") ) {
                speciesTypeKey = SpeciesType.parse(cols[5].substring(8));
            }

            // from-gene validation
            Gene fromGene = null;

            if(false) {
                int fromRgdId = 0;
                List<Gene> genes = dao.getGenesBySymbolAndSpecies(fromGeneSymbol, speciesTypeKey);
                // remove inactive genes
                Iterator<Gene> it = genes.iterator();
                while (it.hasNext()) {
                    Gene gene = it.next();
                    if (!dao.getRgdId(gene.getRgdId()).getObjectStatus().equals("ACTIVE")) {
                        it.remove();
                    }
                }
                // only one active gene must be left
                if (genes.size() != 1) {
                    System.out.println("No valid genes for symbol " + fromGeneSymbol);
                    out.write(line);
                    out.newLine();
                    continue;
                }
                fromRgdId = genes.get(0).getRgdId();

                // the from gene rgd id must match the supplied from gene rgd id
                if (fromRgdIdIncoming != 0) {
                    if (fromRgdId != fromRgdIdIncoming) {
                        System.out.println("Incoming gene rgd id mismatch for " + fromGeneSymbol);
                        out.write(line);
                        out.newLine();
                        continue;
                    }
                }
                // match incoming gene by eg-id
                genes = dao.getGenesByEGID(fromGeneId);
                if( genes.size()!=1 ) {
                    System.out.println("multiple or none genes for from-EGID "+fromGeneSymbol);
                    out.write(line); out.newLine();
                    continue;
                }
                fromGene = genes.get(0);
                if( fromGene.getRgdId()!=fromRgdId ) {
                    System.out.println("EGID1 mismatch for "+fromGeneSymbol);
                    out.write(line); out.newLine();
                    continue;
                }
            } else {
                // match incoming gene by eg-id
                List<Gene> genes = dao.getGenesByEGID(fromGeneId);
                if( genes.size()!=1 ) {
                    System.out.println("multiple or none genes for from-EGID "+fromGeneSymbol);
                    out.write(line); out.newLine();
                    continue;
                }
                fromGene = genes.get(0);
                if (!dao.getRgdId(fromGene.getRgdId()).getObjectStatus().equals("ACTIVE")) {
                    System.out.println("old gene is not active "+fromGeneSymbol);
                    continue;
                }
            }


            // match incoming gene by eg-id
            List<Gene> genes = dao.getGenesByEGID(toGeneId);
            if( genes.size()!=1 ) {
                System.out.println("multiple or none genes for to-EGID "+fromGeneSymbol);
                out.write(line); out.newLine();
                continue;
            }
            Gene toGene = genes.get(0);
            if (!dao.getRgdId(toGene.getRgdId()).getObjectStatus().equals("ACTIVE")) {
                System.out.println("new gene is not active "+toGene.getSymbol());
                continue;
            }

            // to-gene must have from-gene symbol as an alias
            if(false) {
                boolean symbolMatchesAlias = false;
                for (Alias a : dao.getAliases(toGene.getRgdId())) {
                    if (Utils.stringsAreEqualIgnoreCase(a.getValue(), fromGeneSymbol)) {
                        symbolMatchesAlias = true;
                        break;
                    }
                }
                if (!symbolMatchesAlias) {
                    System.out.println("from gene symbol does not match alias in to-gene " + fromGeneSymbol);
                    out.write(line);
                    out.newLine();
                    continue;
                }
            }

            // we have a match!
            mergeMap.put(fromGene.getRgdId(), toGene.getRgdId());
        }
        reader.close();
        out.close();

        System.out.println("to process: "+mergeMap.size());

        String fname = "/tmp/merge_rgd_ids.txt";
        out = new BufferedWriter(new FileWriter(fname));
        out.write("merge-from-rgd-id\tmerge-to-rgd-id");
        out.newLine();
        for( Map.Entry<Integer,Integer> entry: mergeMap.entrySet() ) {
            out.write(entry.getKey()+"\t"+entry.getValue());
            out.newLine();
        }
        out.close();
        simpleMerge(fname, speciesTypeKey);
    }

    /**
     * Created by IntelliJ IDEA.
     * User: mtutaj
     * Date: 7/27/15
     * Time: 10:10 AM
     * <p>
     * reads a list of genes to be merged from an external file, and then runs a merge
     * input file format:
     * merge-from-rgd-id, merge-to-rgd-id, keep-this-ncbi-gene-id
     */
    public static void simpleMerge(String fname, int speciesTypeKey) throws Exception {

        Counters counters = new Counters();

        EGDAO dao = EGDAO.getInstance();

        // load the input file
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = reader.readLine(); // skip header line
        while( (line=reader.readLine())!=null ) {
            String[] cols = line.split("[\\t]", -1);
            int mergeFromRgdId = Integer.parseInt(cols[0]);
            int mergeToRgdId = Integer.parseInt(cols[1]);
            String keepThisGeneId = cols.length>2 ? cols[2].trim() : null;

            // both rgd ids must be active
            RgdId idFrom = dao.getRgdId(mergeFromRgdId);
            if( !idFrom.getObjectStatus().equals("ACTIVE") || idFrom.getSpeciesTypeKey()!= speciesTypeKey ) {
                System.out.println(mergeFromRgdId+" is NOT an active gene or is wrong species!");
                counters.increment("genes skipped");
                continue;
            }
            RgdId idTo = dao.getRgdId(mergeToRgdId);
            if( !idTo.getObjectStatus().equals("ACTIVE") || idTo.getSpeciesTypeKey()!= speciesTypeKey ) {
                System.out.println(mergeToRgdId+" is NOT an active gene or is wrong species!");
                counters.increment("genes skipped");
                continue;
            }

            Gene geneFrom = dao.getGene(mergeFromRgdId);
            Gene geneTo = dao.getGene(mergeToRgdId);

            simpleMerge(geneFrom, geneTo, keepThisGeneId, counters);
            System.out.println();
        }
        reader.close();

        counters.dump();
    }

    public static void simpleMerge(Gene geneFrom, Gene geneTo, String keepThisGeneId, Counters counters) throws Exception {

        EGDAO dao = EGDAO.getInstance();

        int mergeFromRgdId = geneFrom.getRgdId();
        int mergeToRgdId = geneTo.getRgdId();

        System.out.println((counters.get("GENES PROCESSED")+1)
                +".  merge from "+geneFrom.getSymbol()+" RGD:"+mergeFromRgdId
                +" to "+geneTo.getSymbol()+" RGD:"+mergeToRgdId);

        counters.increment("aliases inserted", handleAliases(geneFrom, geneTo, dao));
        counters.increment("notes inserted", handleNotes(mergeFromRgdId, mergeToRgdId));
        counters.increment("references inserted", handleReferences(mergeFromRgdId, mergeToRgdId));
        counters.increment("annots inserted", handleAnnots(mergeFromRgdId, mergeToRgdId));
        counters.increment("xdb ids inserted", handleXdbIds(mergeFromRgdId, mergeToRgdId, dao, keepThisGeneId));
        counters.increment("nomen inserted", handleNomen(geneFrom, geneTo, dao));
        counters.increment("map data inserted", handleMapData(mergeFromRgdId, mergeToRgdId, dao));
        handleHistory(mergeFromRgdId, mergeToRgdId);
        counters.increment("GENES PROCESSED");
    }

    static int handleAliases(Gene geneFrom, Gene geneTo, EGDAO dao) throws Exception {

        // shall from-gene symbol be added as alias to to-gene?
        List<Alias> aliasesFrom = dao.getAliases(geneFrom.getRgdId());
        List<Alias> aliasesTo = dao.getAliases(geneTo.getRgdId());

        Alias alias = new Alias();
        alias.setTypeName("old_gene_symbol");
        alias.setValue(geneFrom.getSymbol());
        alias.setNotes("created by GeneMerge tool on "+new Date());
        aliasesFrom.add(alias);

        // shall from-gene name be added as alias to to-gene?
        alias = new Alias();
        alias.setTypeName("old_gene_name");
        alias.setValue(geneFrom.getName());
        alias.setNotes("created by GeneMerge tool on "+new Date());
        aliasesFrom.add(alias);

        // remove duplicates
        Iterator<Alias> it = aliasesFrom.iterator();
        while( it.hasNext() ) {
            alias = it.next();
            // alias must not be the same as gene name or gene symbol
            if( Utils.stringsAreEqualIgnoreCase(geneTo.getName(), alias.getValue()) ||
                Utils.stringsAreEqualIgnoreCase(geneTo.getSymbol(), alias.getValue())) {
                it.remove();
                continue;
            }

            // weed out duplicate aliases
            for( Alias alias2: aliasesTo ) {
                if( compareTo(alias, alias2)==0 ) {
                    it.remove();
                    break;
                }
            }
        }

        // are there any aliases to insert?
        if( !aliasesFrom.isEmpty() ) {
            for( Alias a: aliasesFrom ) {
                a.setRgdId(geneTo.getRgdId());
            }
            dao.insertAliases(aliasesFrom);

            System.out.println("  aliases inserted: "+aliasesFrom.size());
        }
        return aliasesFrom.size();
    }

    static int compareTo(Alias a1, Alias a2) {

        int r = Utils.stringsCompareToIgnoreCase(a1.getTypeName(), a2.getTypeName());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareToIgnoreCase(a1.getValue(), a2.getValue());
        return r;
    }

    static int handleNotes(int mergeFromRgdId, int mergeToRgdId) throws Exception {

        NotesDAO notesDAO = new NotesDAO();
        List<Note> notesFrom = notesDAO.getNotes(mergeFromRgdId);
        List<Note> notesTo = notesDAO.getNotes(mergeToRgdId);
        int notesInserted = 0;

        // add unique 'from' notes
        for( Note note: notesFrom ) {
            boolean isDuplicated = false;
            for( Note note2: notesTo ) {
                // notes are the same if they have the same type and text
                if( Utils.stringsAreEqualIgnoreCase(note.getNotes(), note2.getNotes()) &&
                    Utils.stringsAreEqualIgnoreCase(note.getNotesTypeName(), note2.getNotesTypeName()) ) {
                    isDuplicated = true;
                    break;
                }
            }
            // insert note if unique
            // to make the note insertable, its key must be made 0, and its rgd id set to target rgd id
            if( !isDuplicated ) {

                note.setRgdId(mergeToRgdId);
                note.setKey(0);
                notesDAO.updateNote(note);
                notesInserted++;
            }
        }

        if( notesInserted>0 )
            System.out.println("  notes inserted: "+notesInserted);
        return notesInserted;
    }

    static int handleReferences(int mergeFromRgdId, int mergeToRgdId) throws Exception {

        AssociationDAO associationDAO = new AssociationDAO();
        List<Reference> refsFrom = associationDAO.getReferenceAssociations(mergeFromRgdId);
        List<Reference> refsTo = associationDAO.getReferenceAssociations(mergeToRgdId);
        int refsInserted = 0;

        for( Reference ref: refsFrom ) {
            boolean isDuplicated = false;
            for( Reference ref2: refsTo ) {
                // curated references are the same if they have the same rgd id
                if( ref.getRgdId()==ref2.getRgdId() ) {
                    isDuplicated = true;
                    break;
                }
            }
            // insert ref if unique
            if( !isDuplicated ) {
                associationDAO.insertReferenceeAssociation(mergeToRgdId, ref.getRgdId());
                refsInserted++;
            }
        }

        if( refsInserted>0 )
            System.out.println("  refs inserted: "+refsInserted);
        return refsInserted;
    }

    static int handleAnnots(int mergeFromRgdId, int mergeToRgdId) throws Exception {

        AnnotationDAO annotationDAO = new AnnotationDAO();
        List<Annotation> annotsFrom = annotationDAO.getAnnotations(mergeFromRgdId);
        List<Annotation> annotsTo = annotationDAO.getAnnotations(mergeToRgdId);
        int annotsInserted = 0;

        for( Annotation from: annotsFrom ) {
            boolean isDuplicated = false;
            for( Annotation to: annotsTo ) {
                if( compareTo(from, to)==0 ) {
                    isDuplicated = true;
                    break;
                }
            }
            // insert annot if unique
            if( !isDuplicated ) {
                from.setAnnotatedObjectRgdId(mergeToRgdId);
                annotationDAO.insertAnnotation(from);
                annotsInserted++;
            }
        }

        if( annotsInserted>0 )
            System.out.println("  annots inserted: "+annotsInserted);
        return annotsInserted;
    }

    // compare by unique key, except annotation object rgd id
    static int compareTo(Annotation a1, Annotation a2) {

        int r = Utils.intsCompareTo(a1.getRefRgdId(), a2.getRefRgdId());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(a1.getTermAcc(), a2.getTermAcc());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(a1.getXrefSource(), a2.getXrefSource());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(a1.getQualifier(), a2.getQualifier());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(a1.getWithInfo(), a2.getWithInfo());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(a1.getEvidence(), a2.getEvidence());
        return r;
    }

    static int handleXdbIds(int mergeFromRgdId, int mergeToRgdId, EGDAO dao, String toBeKeptGeneId) throws Exception {

        XdbId filter = new XdbId();
        filter.setRgdId(mergeFromRgdId);
        List<XdbId> xdbIdsFrom = dao.getXdbIds(filter);
        filter.setRgdId(mergeToRgdId);
        List<XdbId> xdbIdsTo = dao.getXdbIds(filter);

        // if to be kept gene id is not give, use the gene id of merge-to gene
        if( toBeKeptGeneId==null ) {
            for( XdbId id: xdbIdsTo ) {
                if( id.getXdbKey()==XdbId.XDB_KEY_ENTREZGENE ) {
                    toBeKeptGeneId = id.getAccId();
                    break;
                }
            }
        }

        // delete entrez gene ids in 'to' gene that are different than to-be-kept-gene-id
        for( XdbId id: xdbIdsTo ) {
            if( id.getXdbKey()==XdbId.XDB_KEY_ENTREZGENE && !id.getAccId().equals(toBeKeptGeneId) ) {
                ArrayList<XdbId> toBeDeletedXdbIds = new ArrayList<>();
                toBeDeletedXdbIds.add(id);
                dao.deleteXdbIds(toBeDeletedXdbIds, "GENE");
            }
        }

        // delete entrez gene ids in 'from' gene that are different than to-be-kept-gene-id
        for( int i=xdbIdsFrom.size()-1; i>=0; i-- ) {
            XdbId id = xdbIdsFrom.get(i);
            if( id.getXdbKey()==XdbId.XDB_KEY_ENTREZGENE && !id.getAccId().equals(toBeKeptGeneId) ) {
                xdbIdsFrom.remove(i);
            }
        }

        // add unique xdb ids
        List<XdbId> forInsert = new ArrayList<>();
        for( XdbId from: xdbIdsFrom ) {
            boolean isDuplicated = false;
            for( XdbId to: xdbIdsTo ) {
                if( compareTo(from, to)==0 ) {
                    isDuplicated = true;
                    break;
                }
            }
            // insert annot if unique
            if( !isDuplicated ) {
                forInsert.add(from);
            }
        }

        if( !forInsert.isEmpty() ) {
            for( XdbId id: forInsert ) {
                String notes = "created by GeneMerge from RGD ID "+mergeFromRgdId;
                if( id.getNotes()==null )
                    id.setNotes(notes);
                else if( !id.getNotes().contains(notes) ) {
                    id.setNotes(id.getNotes()+"; "+notes);
                }

                id.setRgdId(mergeToRgdId);
                id.setModificationDate(new Date());
            }

            dao.insertXdbs(forInsert, "GENE");
            System.out.println("  xrefs inserted: "+forInsert.size());
        }
        return forInsert.size();
    }

    // compare by unique key, except rgd id
    static int compareTo(XdbId x1, XdbId x2) {

        int r = x1.getXdbKey() - x2.getXdbKey();
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(x1.getAccId(), x2.getAccId());
        return r;
    }

    static int handleNomen(Gene geneFrom, Gene geneTo, EGDAO dao) throws Exception {

        List<NomenclatureEvent> nomenFrom = dao.getNomenclatureEvents(geneFrom.getRgdId());
        List<NomenclatureEvent> nomenTo = dao.getNomenclatureEvents(geneTo.getRgdId());
        int nomenInserted = 0;

        // create nomenclature event for 'to-gene'
        NomenclatureEvent ev = new NomenclatureEvent();
        ev.setDesc("Data Merged");
        ev.setEventDate(new Date());
        ev.setName(geneTo.getName());
        ev.setNomenStatusType("PROVISIONAL");
        ev.setNotes("GeneMerge from RGD ID " + geneFrom.getRgdId() + " to RGD ID " + geneTo.getRgdId());
        ev.setOriginalRGDId(geneFrom.getRgdId());
        ev.setPreviousName(geneFrom.getName());
        ev.setPreviousSymbol(geneFrom.getSymbol());
        ev.setRefKey("9333");
        ev.setRgdId(geneTo.getRgdId());
        ev.setSymbol(geneTo.getSymbol());
        nomenFrom.add(ev);

        // add unique nomen
        for( NomenclatureEvent from: nomenFrom ) {
            boolean isDuplicated = false;
            for( NomenclatureEvent to: nomenTo ) {
                if( compareTo(from, to)==0 ) {
                    isDuplicated = true;
                    break;
                }
            }
            // insert nomen if unique
            if( !isDuplicated ) {
                String notes = "created by GeneMerge from RGD ID "+geneFrom.getRgdId();
                if( from.getNotes()==null )
                    from.setNotes(notes);
                else
                    from.setNotes(from.getNotes()+"; "+notes);

                from.setRgdId(geneTo.getRgdId());

                dao.createNomenEvent(from);
                nomenInserted++;
            }
        }

        if( nomenInserted>0 )
            System.out.println("  nomen inserted: "+nomenInserted);
        return nomenInserted;
    }

    static int compareTo(NomenclatureEvent e1, NomenclatureEvent e2) {

        int r = Utils.stringsCompareTo(e1.getNomenStatusType(), e2.getNomenStatusType());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(e1.getDesc(), e2.getDesc());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(e1.getSymbol(), e2.getSymbol());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(e1.getName(), e2.getName());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(e1.getPreviousSymbol(), e2.getPreviousSymbol());
        if( r!=0 )
            return r;
        r = Utils.stringsCompareTo(e1.getPreviousName(), e2.getPreviousName());
        if( r!=0 )
            return r;
        // compare year, month and day of event date
        r = compareDateTo(e1.getEventDate(), e2.getEventDate());
        return r;
    }

    static int compareDateTo(Date d1, Date d2) {
        int r = d1.getYear() - d2.getYear();
        if( r!=0 )
            return r;
        r = d1.getMonth() - d2.getMonth();
        if( r!=0 )
            return r;
        r = d1.getDay() - d2.getDay();
        return r;
    }

    static void handleHistory(int rgdIdFrom, int rgdIdTo) throws Exception {

        RGDManagementDAO dao = new RGDManagementDAO();
        int rgdId = dao.getRgdIdFromHistory(rgdIdFrom);
        if( rgdId!=rgdIdTo ) {
            dao.recordIdHistory(rgdIdFrom, rgdIdTo);
            System.out.println("  id history: " + rgdIdFrom+" --> "+ rgdIdTo);
        }
        dao.retire(dao.getRgdId2(rgdIdFrom));

        dao.updateLastModifiedDate(rgdIdFrom);
        dao.updateLastModifiedDate(rgdIdTo);
    }

    static int handleMapData(int fromRgdId, int toRgdId, EGDAO dao) throws Exception {

        // if 'to-gene' does not have a positions on specific assemblies of 'from-gene',
        //  then these positions are moved to 'to-gene'
        List<MapData> fromMapData = dao.getMapData(fromRgdId);
        List<MapData> toMapData = dao.getMapData(toRgdId);
        List<MapData> forInsertMapData = new ArrayList<>();

        for( MapData from: fromMapData ) {
            // from-map-key must not be among to-map-keys
            boolean isDuplicate = false;
            for( MapData to: toMapData ) {
                if( to.getMapKey().equals(from.getMapKey()) ) {
                    isDuplicate = true;
                    break;
                }
            }
            if( !isDuplicate ) {
                forInsertMapData.add(from);
            }
        }
        if( !forInsertMapData.isEmpty() ) {
            for( MapData md: forInsertMapData ) {
                md.setRgdId(toRgdId);
            }
            dao.insertMapData(forInsertMapData);
            System.out.println("  map data inserted: " + forInsertMapData.size());
        }
        return forInsertMapData.size();
    }
}
