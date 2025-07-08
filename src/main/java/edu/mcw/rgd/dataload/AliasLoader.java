package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.EGDAO;
import edu.mcw.rgd.datamodel.Alias;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;

/**
 * @author mtutaj
 * @since 3/25/14
 * code for handling aliases (synonyms)
 */
public class AliasLoader {

    protected final Logger logAliases = LogManager.getLogger("aliases");

    private List<Alias> incoming = new ArrayList<>(); // aliases retrieved from eg gene record
    public List<Alias> forInsert; // new aliases to be added to rgd db (subset of 'aliases')

    public List<Alias> getIncoming() {
        return incoming;
    }

    public void addAlias(Alias alias) {

        // never add an alias without a value
        if( Utils.isStringEmpty(alias.getValue()) ) {
            return;
        }

        // do not add 'null','test','unknown' aliases
        if( Utils.stringsAreEqualIgnoreCase("null", alias.getValue()) ||
            Utils.stringsAreEqualIgnoreCase("test", alias.getValue()) ||
            Utils.stringsAreEqualIgnoreCase("unknown", alias.getValue()) )
            return;

        // strip potential prefix: 'LOW QUALITY PROTEIN:'
        String prefix = "LOW QUALITY PROTEIN:";
        if( alias.getValue().startsWith(prefix) ) {
            alias.setValue( alias.getValue().substring(prefix.length()).trim() );
        }

        // compress whitespace in value: replace multiple space characters with a single space
        compressWhitespace(alias);

        // make sure the to-be-added alias is unique
        for( Alias a: incoming ) {
            if( Utils.stringsAreEqualIgnoreCase(a.getValue(), alias.getValue()) ) {
                // alias 'alias' is already in the array
                return;
            }
        }

        if( !handleSemicolons(alias) )
            incoming.add(alias);
    }

    void compressWhitespace(Alias alias) {
        alias.setValue(alias.getValue().replaceAll("[\\s]+", " "));
    }

    boolean handleSemicolons(Alias a) {

        // no special processing if no semicolons
        if( !a.getValue().contains("; ") )
            return false;

        // examine every part: its parentheses '(',')' must match
        String[] parts = a.getValue().split("; ");
        for( String part: parts ) {
            int parentheseCount = 0;
            for( int i=0; i<part.length(); i++ ) {
                char c = part.charAt(i);
                if( c=='(' )
                    parentheseCount++;
                else if( c==')' )
                    parentheseCount--;
            }
            if( parentheseCount!=0 ) {
                // problems with parentheses, f.e. "(NAD+; poly (ADP-ribose) polymerase)"
                // these conflicts are ported by qc
                //   1 solution is that a manual equivalent synonym is created,
                //                            f.e. "(NAD+, poly (ADP-ribose) polymerase)"
                String equivalentValue = a.getValue().replace("; ", ", ");
                a.setValue(equivalentValue);
                addAlias(a);
                return false;
            }
        }

        // parentheses are OK: create alias for every part
        for( String part: parts ) {
            Alias a2 = new Alias();
            a2.setNotes(a.getNotes());
            a2.setRgdId(a.getRgdId());
            a2.setSpeciesTypeKey(a.getSpeciesTypeKey());
            a2.setTypeName(a.getTypeName());
            a2.setValue(part);
            addAlias(a2);
        }
        return true;
    }

    // qc aliases if gene is to be updated
    public void qcIncomingAliases(int rgdId, String symbol, String name, CounterPool counters) throws Exception {

        // incoming aliases: removed protein-redundant aliases
        removeProteinRedundantAliases(incoming, symbol, rgdId);
        removeProteinRedundantAliases(incoming, name, rgdId);

        // update alias, keep the old alias and add new alias
        List<Alias> rgdAliases = EGDAO.getInstance().getAliases(rgdId);
        forInsert = new ArrayList<>(); // new aliases to be added to ALIASES table

        for (Alias alias: incoming) {
            if( isAliasOnList(rgdAliases, alias)) {
                //logger.debug("gene alias "+ alias.getValue() +" is already in RGD");
                counters.increment("ALIASES_MATCHED");
            }
            else {
                alias.setRgdId(rgdId);
                forInsert.add(alias);
                logAliases.debug("gene alias ["+ alias.getValue() +"] inserted into alias table");
            }
        }
    }

    // remove those aliases from the to-be-inserted aliases that happen to be the same as gene name or symbol
    public void qcAliases(BulkGene bg, CounterPool counters) {

        if( forInsert==null || forInsert.isEmpty() )
            return;

        // create a set of taboo aliases
        Set<String> tabooAliasesLC = new HashSet<>();
        Gene gene = bg.getGene();
        int rgdId = 0;
        if( gene!=null ) {
            rgdId = gene.getRgdId();
            if( gene.getSymbol()!=null )
                tabooAliasesLC.add(gene.getSymbol().toLowerCase());
            if( gene.getName()!=null )
                tabooAliasesLC.add(gene.getName().toLowerCase());
        }
        gene = bg.rgdGene;
        if( gene!=null ) {
            if( rgdId==0 )
                rgdId = gene.getRgdId();
            if( gene.getSymbol()!=null )
                tabooAliasesLC.add(gene.getSymbol().toLowerCase());
            if( gene.getName()!=null )
                tabooAliasesLC.add(gene.getName().toLowerCase());
        }

        // remove all taboo aliases
        Iterator<Alias> it = forInsert.iterator();
        while( it.hasNext() ) {
            String aliasValue = it.next().getValue();
            if( aliasValue!=null && tabooAliasesLC.contains(aliasValue.toLowerCase()) ) {
                // skip this alias if nomenclature event is not in place
                // (else we will skip alias that is tracking gene name/symbol change!)
                if( !bg.isFlagSet("ORTHO_NOMEN_NAME") && !bg.isFlagSet("ORTHO_NOMEN_SYMBOL") ){
                    counters.increment("SKIPPED_ALIASES_SAME_AS_GENE_NAME_SYMBOL");
                    it.remove();
                    logAliases.debug("*** RGDID:" + rgdId + " skipped alias, same as gene symbol/name[" + aliasValue + "]");
                }
            }
        }
    }

    // return true if given alias value is on the list
    boolean isAliasOnList(List<Alias> aliases, Alias newAlias) {
        for( Alias alias:aliases ) {
            if( Utils.stringsAreEqualIgnoreCase(alias.getValue(), newAlias.getValue()))
                return true;
        }
        return false;
    }

    /**
     * incoming aliases cannot have redundant entries regarding gene symbol and 'protein' suffix/prefix;
     * f.e. if gene symbol is 'Fam216a', aliases 'protein Fam216a' and 'Fam216 protein' should be removed;
     *  similarly if there are aliases 'Rbs11' and 'Rbs11 protein', 'Rbs11 protein' should be removed as well;
     *  similarly if there are aliases 'Rbs11' and 'protein Rbs11', 'protein Rbs11' should be removed as well
     * @param aliases list of incoming aliases
     * @param geneSymbol gene symbol
     * @param rgdId rgd id
     */
    void removeProteinRedundantAliases(List<Alias> aliases, String geneSymbol, int rgdId) {

        String lcGeneSymbol = Utils.defaultString(geneSymbol).toLowerCase();
        if( lcGeneSymbol.isEmpty() ) {
            return;
        }

        // if alias is prefixed/suffixed by 'protein', the remaining part is held here
        Set<String> possibleDuplicates = new HashSet<>();

        // remove duplicate symbols
        Iterator<Alias> it = aliases.iterator();
        while( it.hasNext() ) {
            Alias a = it.next();
            String lcAlias = a.getValue().toLowerCase();
            if( lcAlias.startsWith("protein ") ) {
                String possibleDuplicate = lcAlias.substring(8);
                possibleDuplicates.add(possibleDuplicate);
                if( possibleDuplicate.equals(lcGeneSymbol) ) {
                    logAliases.debug("*** RGDID:" + rgdId + " skipped alias duplicate [" + a.getValue() + "]");
                    it.remove();
                }
                continue;
            }
            if( lcAlias.endsWith(" protein") ) {
                String possibleDuplicate = lcAlias.substring(0, lcAlias.length()-8);
                possibleDuplicates.add(possibleDuplicate);
                if( possibleDuplicate.equals(lcGeneSymbol) ) {
                    logAliases.debug("*** RGDID:" + rgdId + " skipped alias duplicate [" + a.getValue() + "]");
                    it.remove();
                }
            }
        }

        // if there are any possible duplicates, then check all aliases against them
        if( !possibleDuplicates.isEmpty() ) {

            for( String possibleDuplicateStem: possibleDuplicates ) {
                it = aliases.iterator();
                while( it.hasNext() ) {
                    Alias a = it.next();
                    if( Utils.stringsAreEqualIgnoreCase(a.getValue(), possibleDuplicateStem+" protein") ||
                        Utils.stringsAreEqualIgnoreCase(a.getValue(), "protein "+possibleDuplicateStem))
                    {
                        logAliases.debug("*** RGDID:" + rgdId + " skipped alias duplicate [" + a.getValue() + "]");
                        it.remove();
                    }
                }
            }
        }
    }

}
