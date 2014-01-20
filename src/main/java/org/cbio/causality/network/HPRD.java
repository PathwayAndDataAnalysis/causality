package org.cbio.causality.network;

import org.cbio.causality.analysis.Graph;

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
		graph = new Graph();

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream(
			"BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[3])) continue;

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
		Graph graph = HPRD.getGraph(true);
		System.out.println("graph.getNeighbors(\"TP53\").size() = " +
			graph.getNeighbors("TP53").size());
	}
}
