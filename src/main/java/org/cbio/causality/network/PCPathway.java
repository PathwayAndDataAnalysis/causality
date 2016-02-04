package org.cbio.causality.network;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.EntityReference;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.Xref;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * @author Ozgun Babur
 */
public class PCPathway
{
	private static Map<String, Set<String>> gene2pathway;
	private static Map<String, Set<String>> pathway2gene;
	private static Map<String, String> pathway2name;

	private static final String FILE = "pcpathway.txt";

	public static Set<String> getPathways(String gene)
	{
		if (gene2pathway.containsKey(gene)) return gene2pathway.get(gene);
		else return Collections.emptySet();
	}

	public static Set<String> getGenes(String pathwayID)
	{
		if (pathway2gene.containsKey(pathwayID)) return pathway2gene.get(pathwayID);
		else return Collections.emptySet();
	}

	public static String getName(String id)
	{
		return pathway2name.get(id);
	}

	public static String getCoverageStr(String id, Set<String> genes)
	{
		if (!pathway2gene.containsKey(id)) return null;
		Set<String> set = new HashSet<String>(genes);
		set.retainAll(pathway2gene.get(id));
		return set.size() + "/" + pathway2gene.get(id).size();
	}

	/**
	 * Gets the enrichment pvals and and pval limits.
	 * @param genes query
	 * @param background if there is any
	 * @return two maps, first is for pvals, second is for limits
	 */
	public static Map<String, Double>[] getEnrichmentPvals(Collection<String> genes,
		Collection<String> background, int minMemberSize)
	{
		if (background == null)
		{
			background = new HashSet<String>(gene2pathway.keySet());
			if (!background.containsAll(genes))
			{
				Set<String> set = new HashSet<String>(genes);
				set.removeAll(background);

				genes = new HashSet<String>(genes);
				genes.removeAll(set);
				System.out.println("Removed " + set.size() + " unknown genes: " + set);
				System.out.println("Using " + genes.size() + ": " + genes);
			}
		}

		if (!background.containsAll(genes)) throw new IllegalArgumentException(
			"Background genes have to contain all the selected genes.");

		Map<String, Integer> selectionCnt = count(genes);
		Map<String, Integer> backgroundCnt = count(background);

		Map<String, Double> mapP = new HashMap<String, Double>();
		Map<String, Double> mapL = new HashMap<String, Double>();

		for (String pathway : selectionCnt.keySet())
		{
			if (pathway2gene.get(pathway).size() < minMemberSize) continue;

			double pval = FishersExactTest.calcEnrichmentPval(background.size(),
				backgroundCnt.get(pathway), genes.size(), selectionCnt.get(pathway));

			double limit = FishersExactTest.calcEnrichmentPval(background.size(),
				backgroundCnt.get(pathway), genes.size(),
				Math.min(backgroundCnt.get(pathway), genes.size()));

			mapP.put(pathway, pval);
			mapL.put(pathway, limit);
		}

		return new Map[]{mapP, mapL};
	}

	private static Map<String, Integer> count(Collection<String> genes)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();

		for (String pathway : pathway2gene.keySet())
		{
			Set<String> mems = new HashSet<String>(pathway2gene.get(pathway));
			mems.retainAll(genes);
			if (!mems.isEmpty()) cnt.put(pathway, mems.size());
		}
		return cnt;
	}

	public static List<String> getEnrichedPathways(Collection<String> genes,
		Collection<String> background, double fdrThr)
	{
		return getEnrichedPathways(genes, background, fdrThr, 3);
	}

	public static List<String> getEnrichedPathways(Collection<String> genes,
		Collection<String> background, double fdrThr, int memberThreshold)
	{
		Map<String, Double>[] map = getEnrichmentPvals(genes, background, memberThreshold);
		if (fdrThr < 0)
		{
			fdrThr = FDR.decideBestFDR_BH(map[0], map[1]);
			System.out.println("fdrThr = " + fdrThr);
		}
		return FDR.select(map[0], map[1], fdrThr);
	}

	/**
	 * Gets pathways sorted to their containment of the query genes. Does not control the density,
	 * but only the count, so it is likely that first pathways will be big ones.
	 */
	public static TreeMap<String, Integer> getSortedPathways(Collection<String> genes)
	{
		final Map<String, Integer> map = new HashMap<String, Integer>();

		for (String gene : genes)
		{
			for (String pathway : getPathways(gene))
			{
				if (map.containsKey(pathway)) map.put(pathway, map.get(pathway) + 1);
				else map.put(pathway, 1);
			}
		}

		TreeMap<String, Integer> sorted = new TreeMap<String, Integer>(new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return map.get(o2).compareTo(map.get(o1));
			}
		});

		sorted.putAll(map);
		return sorted;
	}

	static
	{
		readResources();
	}

	public static void main(String[] args) throws IOException
	{
//		String s = "PTEN, PIK3CA, ARID1A, PIK3R1, CTNNB1, TP53, KRAS, CTCF, FBXW7, LRP2, FGFR2, RYR1, TBL1XR1, MTOR, CACNA1A, PPP2R1A, PKN1, LYST, TRPM6, ERBB2, FN1, WDFY3, MYC, SPTB, DVL3, PRKCI, ECT2, ACTL6A, TBC1D31, IKBKB, PRKACA, DLG1, PTK2, THPO, DNM2, FOSL2, DSTYK, CCNE1, TNK2, EFNA1, PAK2, RASAL1, ARMC6, HGS, CDC37, TNFSF10, PPP1R1B, GRB2, PPP1CA";
//		List<String> select = getEnrichedPathways(Arrays.asList(s.split(", ")), null, 0.01);
//		for (String id : select)
//		{
//			System.out.println(id + "\t" + getName(id));
//		}

		prepareResource();
	}

	private static void prepareResource() throws IOException
	{
		pathway2gene = new HashMap<String, Set<String>>();
		gene2pathway = new HashMap<String, Set<String>>();
		pathway2name = new HashMap<String, String>();

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("Pathway Commons.7.Detailed_Process_Data.BIOPAX.owl"));
		for (Pathway pathway : model.getObjects(Pathway.class))
		{
			String id = pathway.getRDFId();
			String name = pathway.getDisplayName();

			if (name == null || name.isEmpty()) continue;

			pathway2name.put(id, name);

			Model m = excise(model, pathway);
			Set<String> syms = collectGeneSymbols(m);

			if (syms.isEmpty()) continue;

			if (!pathway2gene.containsKey(id)) pathway2gene.put(id, new HashSet<String>());

			for (String sym : syms)
			{
				if (!gene2pathway.containsKey(sym)) gene2pathway.put(sym, new HashSet<String>());

				pathway2gene.get(id).add(sym);
				gene2pathway.get(sym).add(id);
			}
		}
		writeResources();
	}

	private static void writeResources() throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("src/main/resources/org/cbio/causality/network/" + FILE));

		for (String id : pathway2gene.keySet())
		{
			if (pathway2gene.get(id).isEmpty()) continue;

			writer.write(id + "\t" + pathway2name.get(id));

			for (String sym : pathway2gene.get(id))
			{
				writer.write("\t" + sym);
			}
			writer.write("\n");
		}
		writer.close();
	}

	private static void readResources()
	{
		pathway2gene = new HashMap<String, Set<String>>();
		gene2pathway = new HashMap<String, Set<String>>();
		pathway2name = new HashMap<String, String>();

		Scanner sc = new Scanner(PCPathway.class.getResourceAsStream(FILE));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");

			pathway2name.put(token[0], token[1]);

			if (token.length > 2)
			{
				pathway2gene.put(token[0],
					new HashSet<String>(Arrays.asList(token).subList(2, token.length)));
			}
			for (int i = 2; i < token.length; i++)
			{
				if (!gene2pathway.containsKey(token[i])) gene2pathway.put(token[i], new HashSet<String>());
				gene2pathway.get(token[i]).add(token[0]);
			}
		}
	}

	private static Model excise(Model model, Pathway pathway)
	{
		Completer c = new Completer(SimpleEditorMap.L3);

		Set<BioPAXElement> objects = c.complete(Collections.<BioPAXElement>singleton(pathway), model);

		Cloner cln = new Cloner(SimpleEditorMap.L3, BioPAXLevel.L3.getDefaultFactory());

		return cln.clone(model, objects);
	}

	private static Set<String> collectGeneSymbols(Model model)
	{
		Set<String> symbols = new HashSet<String>();

		for (EntityReference er : model.getObjects(EntityReference.class))
		{
			for (Xref xref : er.getXref())
			{
				if (xref.getDb() == null) continue;
				if (xref.getDb().equalsIgnoreCase("HGNC SYMBOL"))
				{
					String s = HGNC.getSymbol(xref.getId());
					if (s != null) symbols.add(s);
				}
			}
		}
		return symbols;
	}
}
