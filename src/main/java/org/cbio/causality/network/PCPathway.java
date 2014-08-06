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
import org.cbio.causality.util.FishersExactTest;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PCPathway
{
	private static Map<String, Set<String>> gene2pathway;
	private static Map<String, Set<String>> pathway2gene;
	private static final String FILE = "pcpathway.txt";

	public static Set<String> getPathways(String gene)
	{
		if (gene2pathway.containsKey(gene)) return gene2pathway.get(gene);
		else return Collections.emptySet();
	}

	public static Set<String> getGenes(String pathwayName)
	{
		if (pathway2gene.containsKey(pathwayName)) return pathway2gene.get(pathwayName);
		else return Collections.emptySet();
	}

	public static Map<String, Double> getEnrichedPathways(Collection<String> genes,
		Collection<String> background)
	{
		if (!background.containsAll(genes)) throw new IllegalArgumentException(
			"Background genes have to contain all the selected genes.");

		Map<String, Integer> selectionCnt = count(genes);
		Map<String, Integer> backgroundCnt = count(background);

		Map<String, Double> map = new HashMap<String, Double>();

		for (String pathway : selectionCnt.keySet())
		{
			double pval = FishersExactTest.calcEnrichmentPval(background.size(),
				backgroundCnt.get(pathway), genes.size(), selectionCnt.get(pathway));

			map.put(pathway, pval);
		}

		return map;
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
		String s = "PTEN, PIK3CA, ARID1A, PIK3R1, CTNNB1, TP53, KRAS, CTCF, FBXW7, LRP2, FGFR2, RYR1, TBL1XR1, MTOR, CACNA1A, PPP2R1A, PKN1, LYST, TRPM6, ERBB2, FN1, WDFY3, MYC, SPTB, DVL3, PRKCI, ECT2, ACTL6A, TBC1D31, IKBKB, PRKACA, DLG1, PTK2, THPO, DNM2, FOSL2, DSTYK, CCNE1, TNK2, EFNA1, PAK2, RASAL1, ARMC6, HGS, CDC37, TNFSF10, PPP1R1B, GRB2, PPP1CA";
		TreeMap<String, Integer> map = getSortedPathways(Arrays.asList(s.split(", ")));
		System.out.println(map);

//		prepareResource();
	}

	private static void prepareResource() throws IOException
	{
		pathway2gene = new HashMap<String, Set<String>>();
		gene2pathway = new HashMap<String, Set<String>>();

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("../biopax-pattern/All-Data.owl"));
		for (Pathway pathway : model.getObjects(Pathway.class))
		{
			String name = pathway.getDisplayName();
			if (name == null || name.isEmpty()) continue;

			Model m = excise(model, pathway);
			Set<String> syms = collectGeneSymbols(m);

			if (syms.isEmpty()) continue;

			if (!pathway2gene.containsKey(name)) pathway2gene.put(name, new HashSet<String>());

			for (String sym : syms)
			{
				if (!gene2pathway.containsKey(sym)) gene2pathway.put(sym, new HashSet<String>());

				pathway2gene.get(name).add(sym);
				gene2pathway.get(sym).add(name);
			}
		}
		writeResources();
	}

	private static void writeResources() throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("src/main/resources/org/cbio/causality/network/" + FILE));

		for (String name : pathway2gene.keySet())
		{
			if (pathway2gene.get(name).isEmpty()) continue;

			writer.write(name);

			for (String sym : pathway2gene.get(name))
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

		Scanner sc = new Scanner(PCPathway.class.getResourceAsStream(FILE));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");

			if (token.length > 1)
			{
				pathway2gene.put(token[0],
					new HashSet<String>(Arrays.asList(token).subList(1, token.length)));
			}
			for (int i = 1; i < token.length; i++)
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
				if (xref.getDb().equals("HGNC Symbol")) symbols.add(xref.getId());
			}
		}
		return symbols;
	}
}
