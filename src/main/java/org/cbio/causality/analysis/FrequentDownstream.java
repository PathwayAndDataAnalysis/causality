package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.idmapping.CancerGeneCensus;
import org.cbio.causality.network.PathwayCommons;

import java.util.*;

/**
 * Given a set of genes, finds the frequent downstream genes.
 * @author Ozgun Babur
 */
public class FrequentDownstream
{
	public static Map<String, Double> getScoredDownstream(Set<String> genes, Graph graph, int depth)
	{
		Map<String, Double> scores = new HashMap<String, Double>();
		for (String gene : genes)
		{
			List<Set<String>> tiers = graph.getNeighborsTiered(Collections.singleton(gene), depth, false);
			for (int i = 0; i < tiers.size(); i++)
			{
				Set<String> tier = tiers.get(i);
				for (String dwn : tier)
				{
					double score = 1D / (double) (i + 1);
					if (scores.containsKey(dwn)) scores.put(dwn, scores.get(dwn) + score);
					else scores.put(dwn, score);
				}
			}
		}

		for (String gene : genes)
		{
			scores.remove(gene);
		}

		return scores;
	}

	public static void main(String[] args)
	{
		Graph graph = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		graph.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));

		Set<String> genes = new HashSet<String>(Arrays.asList(("FYB, HSP90AA1, LYN, BCAR1, ACP1, PXN, PTPN11, MAPK1, SHB, PLCG1, FYN, MAPK14, MAPK3, LCK, PDGFRA, SHC1, YES1, FRS2").split(", ")));

		final Map<String, Double> map = getScoredDownstream(genes, graph, 5);
		List<String> downstream = new ArrayList<String>(map.keySet());
		Collections.sort(downstream, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return map.get(o2).compareTo(map.get(o1));
			}
		});

		Set<String> census = CancerGeneCensus.getAllSymbols();
		for (String dwn : downstream)
		{
			if (census.contains(dwn))
				System.out.println(graph.getUpstream(dwn).size() + "\t" + dwn + "\t" + map.get(dwn));
		}
		System.out.println("downstream: " + downstream);

		genes.retainAll(census);
		System.out.println(new ArrayList<String>(genes));
	}
}
