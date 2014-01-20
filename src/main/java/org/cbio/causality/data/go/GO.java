package org.cbio.causality.data.go;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.util.FishersExactTest;

import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

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
		Graph graph = new Graph();

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

	public static void main(String[] args)
	{
		Graph graph = getGraph(Namespace.biological_process);
		System.out.println(graph.getDownstream("CELL_JUNCTION"));
		System.out.println(graph.getUpstream("TP53"));
	}
}
