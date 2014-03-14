package org.cbio.causality.data.go;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;
import org.cbio.causality.util.FormatUtil;

import java.util.*;

/**
 * For finding enriched GO terms.
 * @author Ozgun Babur
 */
public class GO
{
	private static Map<Namespace, Graph> graphMap = new HashMap<Namespace, Graph>();

	public static Graph getGraph(Namespace ns)
	{
		if (!graphMap.containsKey(ns)) load(ns);
		return graphMap.get(ns);
	}

	public static Set<String> getMembers(String term, Namespace ns)
	{
		if (!graphMap.containsKey(ns)) load(ns);
		return graphMap.get(ns).getDownstream(term);
	}

	public static Set<String> getTerms(String gene, Namespace ns)
	{
		if (!graphMap.containsKey(ns)) load(ns);
		return graphMap.get(ns).getUpstream(gene);
	}

	public static Map<String, Double> getEnrichedTerms(Set<String> selectedGenes,
		Set<String> backgroundGenes, Namespace ns)
	{
		if (!backgroundGenes.containsAll(selectedGenes)) throw new IllegalArgumentException(
			"Background genes have to contain all the selected genes.");

		if (!graphMap.containsKey(ns)) load(ns);

		Map<String, Integer> selectionCnt = count(selectedGenes, ns);
		Map<String, Integer> backgroundCnt = count(backgroundGenes, ns);

		Map<String, Double> map = new HashMap<String, Double>();

		for (String term : selectionCnt.keySet())
		{
			double pval = FishersExactTest.calcEnrichmentPval(backgroundGenes.size(),
				backgroundCnt.get(term), selectedGenes.size(), selectionCnt.get(term));

			map.put(term, pval);
		}

		return map;
	}

	private static Map<String, Integer> count(Set<String> genes, Namespace ns)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();

		for (String gene : genes)
		{
			for (String term : graphMap.get(ns).getUpstream(gene))
			{
				if (!cnt.containsKey(term)) cnt.put(term, 1);
				else cnt.put(term, cnt.get(term) + 1);
			}
		}

		return cnt;
	}

	private static void load(Namespace ns)
	{
		Graph graph = new Graph("GO", ns.name());

		Scanner sc = new Scanner(GO.class.getResourceAsStream(ns.filename));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			for (int i = 2; i < token.length; i++)
			{
				graph.putRelation(token[0], token[i], true);
			}
		}
		
		graphMap.put(ns, graph);
	}

	public enum Namespace
	{
		any("msigdb-go.gmt"),
		molecular_function("msigdb-go-mf.gmt"),
		biological_process("msigdb-go-bp.gmt"),
		cellular_component("msigdb-go-cc.gmt");

		String filename;

		Namespace(String filename)
		{
			this.filename = filename;
		}
	}

	public static void printEnrichment(Set<String> selectedGenes, Set<String> backgroundGenes,
		double fdrThr)
	{
		for (GO.Namespace ns : GO.Namespace.values())
		{
			if (ns == Namespace.any) continue;

			System.out.println("\n---------------------Go terms - " + ns.name());

			Map<String, Double> goMap = GO.getEnrichedTerms(
				selectedGenes, backgroundGenes, ns);

			List<String> enrichedGO = FDR.select(goMap, null, fdrThr);
			for (String go : enrichedGO)
			{
				List<String> members = new ArrayList<String>(GO.getMembers(go, ns));
				members.retainAll(selectedGenes);
				Collections.sort(members);

				System.out.println(go + "\t" +
					FormatUtil.roundToSignificantDigits(goMap.get(go), 2) + "\t" + members);
			}
		}
	}

	public static void main(String[] args)
	{
		Set<String> terms = getTerms("CCNH", Namespace.any);
		for (String term : terms)
		{
			System.out.println(term);
		}
	}
}
