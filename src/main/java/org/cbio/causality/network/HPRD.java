package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class HPRD implements InteractionProvider
{
	private static Graph graph;

	public static Set<String> getInteractors(String symbol)
	{
		return graph.getNeighbors(symbol);
	}

	public static Set<String> getAllSymbols()
	{
		return graph.getSymbols();
	}

	public static int getDegree(String symbol)
	{
		return graph.getNeighbors(symbol).size();
	}

	static
	{
		graph = new Graph("HPRD", SIFEnum.INTERACTS_WITH.getTag());

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream(
			"BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[3])) continue;

			token[0] = HGNC.getSymbol(token[0]);
			token[3] = HGNC.getSymbol(token[3]);

			if (token[0] == null || token[3] == null) continue;

			graph.putRelation(token[0], token[3], false);
		}
	}

	@Override
	public Set<String> getInteractions(String symbol)
	{
		return HPRD.getInteractors(symbol);
	}

	public static Graph getGraph()
	{
		return graph;
	}

	public static Graph getGraph(boolean directed)
	{
		if (!directed) return getGraph();

		Graph graph = new Graph();
		for (String s : getAllSymbols())
		{
			for (String n : getInteractors(s))
			{
				graph.putRelation(s, n, true);
			}
		}

		return graph;
	}

	public static void main(String[] args)
	{
		Graph icw = PathwayCommons.getGraph(SIFEnum.IN_COMPLEX_WITH);
		Graph hprd = HPRD.getGraph();
		Graph intact = IntAct.getGraph();
		icw.printVennIntersections(false, hprd, intact);
	}
}
