package edu.mcw.rgd.process;

import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.MappedGene;

import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * Created by mtutaj on 8/30/2016.
 */
public class BulkRenamerForTRnaGenes {

    public static void main(String[] args) throws Exception {

        int speciesTypeKey = 3;
        int mapKeyRef = 360;
        int mapKeyCelera = 15;

        GeneDAO geneDAO = new GeneDAO();
        List<Gene> newGeneGenes = geneDAO.getActiveGenes("newgene", speciesTypeKey);
        System.out.println("NEWGENE genes "+newGeneGenes.size());

        // limit results to 'trna' genes
        limitGenesToTrna(newGeneGenes);
        System.out.println("NEWGENE trna genes "+newGeneGenes.size());

        List<MappedGene> mappedGenesRef = geneDAO.getActiveMappedGenesByGeneList(mapKeyRef, newGeneGenes);
        System.out.println("NEWGENE trna ref genes "+mappedGenesRef.size());

        List<MappedGene> mappedGenesCelera = geneDAO.getActiveMappedGenesByGeneList(mapKeyCelera, newGeneGenes);
        System.out.println("NEWGENE trna celera genes "+mappedGenesCelera.size());

        combineGenes(mappedGenesRef, mappedGenesCelera);

        Collections.sort(mappedGenesRef, new Comparator<MappedGene>() {
            @Override
            public int compare(MappedGene o1, MappedGene o2) {
                int r = Utils.stringsCompareTo(o1.getGene().getSymbol(), o2.getGene().getSymbol());
                if( r!=0 )
                    return r;
                // REF genes go before CELERA genes
                r = Utils.stringsCompareTo(o2.getGene().getNotes(), o1.getGene().getNotes());
                if( r!=0 )
                    return r;
                r = Chromosome.getOrdinalNumber(o1.getChromosome()) - Chromosome.getOrdinalNumber(o2.getChromosome());
                if( r!=0 )
                    return r;
                long l = o1.getStart() - o2.getStart();
                if( l!=0 )
                    return l<0 ? -1 : l>0 ? 1 : 0;
                l = o1.getStop() - o2.getStop();
                return l<0 ? -1 : l>0 ? 1 : 0;
            }
        });

        processGeneGroup(mappedGenesRef, geneDAO);
    }

    static void limitGenesToTrna(List<Gene> genes) {

        Iterator<Gene> it = genes.iterator();
        while( it.hasNext() ) {
            Gene gene = it.next();
            if( !gene.getType().equals("trna") ) {
                it.remove();
            }
        }
    }

    static void processGeneGroup(List<MappedGene> genes, GeneDAO geneDAO) throws Exception {

        // gene group has the same symbol
        int geneGroup = 0;
        String groupSymbol = "";
        Gene gene = null;
        for( MappedGene mgene: genes ) {
            System.out.println("RGD: "+mgene.getGene().getRgdId()+" "+mgene.getGene().getSymbol()+" chr"+mgene.getChromosome()+":"+mgene.getStart()+".."+mgene.getStop()+"("+mgene.getStrand()+")");
            // detect new gene group
            if( !Utils.stringsAreEqualIgnoreCase(mgene.getGene().getSymbol(), groupSymbol) ) {
                groupSymbol = mgene.getGene().getSymbol();
                gene = geneDAO.getGene(Integer.parseInt(groupSymbol.substring(8)));
                geneGroup = 0;
            }

            // change gene symbol to the source gene symbol
            String newGeneSymbol = gene.getSymbol()+(++geneGroup);
            String newGeneName = mgene.getGene().getName()+" "+geneGroup;
            System.out.println("  SYMBOL: "+mgene.getGene().getSymbol()+" ==> "+newGeneSymbol+" "+mgene.getGene().getNotes());
            System.out.println("    NAME: "+mgene.getGene().getName()+" ==> "+newGeneName);

            // do the update
            mgene.getGene().setSymbol(newGeneSymbol);
            mgene.getGene().setName(newGeneName);
            //geneDAO.updateGene(mgene.getGene());
        }
    }

    static void combineGenes(List<MappedGene> genesRef, List<MappedGene> genesCelera) {
        for( MappedGene rg: genesRef ) {
            rg.getGene().setNotes("REF");
        }

        // add to ref genes, celera only genes
        for( MappedGene cg: genesCelera ) {
            boolean match = false;
            for( MappedGene rg: genesRef ) {
                if( rg.getGene().getRgdId()==cg.getGene().getRgdId() ) {
                    match = true;
                    break;
                }
            }
            if( !match ) {
                // celera gene does not match ref gene
                cg.getGene().setNotes("CELERA");
                genesRef.add(cg);
            }
        }
    }
}
