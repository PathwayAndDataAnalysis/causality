package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.io.*;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class IntAct implements InteractionProvider
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
		graph = new Graph("IntAct", SIFEnum.INTERACTS_WITH.getTag());

		Scanner sc = new Scanner(IntAct.class.getResourceAsStream("IntAct.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[1])) continue;

			token[0] = HGNC.getSymbol(token[0]);
			token[1] = HGNC.getSymbol(token[1]);

			if (token[0] == null || token[1] == null) continue;

			graph.putRelation(token[0], token[1], false);
		}
	}

	@Override
	public Set<String> getInteractions(String symbol)
	{
		return IntAct.getInteractors(symbol);
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

	public static void main(String[] args) throws IOException
	{
		Graph intact = getGraph();
		intact.printVennIntersections(HPRD.getGraph());
		intact.printVennIntersections(PathwayCommons.getGraph(SIFEnum.IN_COMPLEX_WITH));
	}

	//--- Section: Preparing IntAct.txt file ------------------------------------------------------|

	private static void prepareResource() throws IOException
	{
		String dir = "/home/ozgun/Downloads/";

		Scanner sc = new Scanner(new File(dir + "intact.txt"));


		Set<String> all = new HashSet<String>();
		Set<String> half = new HashSet<String>();

		while(sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			String g1 = extractGeneSymbol(token[4]);
			String g2 = extractGeneSymbol(token[5]);

			if (g1 != null && g2 != null)
			{
				line = g1 + "\t" + g2;

				if (all.contains(line)) continue;
				String rev = g2 + "\t" + g1;
				assert !all.contains(rev);

				all.add(line);
				half.add(line);
				all.add(rev);
			}
		}

		System.out.println("size = " + half.size());

		BufferedWriter writer = new BufferedWriter(new FileWriter("IntAct.txt"));

		for (String s : half)
		{
			writer.write(s + "\n");
		}

		writer.close();
	}

	private static String extractGeneSymbol(String s)
	{
		int i = s.indexOf("(gene name)");

		if (i > 0)
		{
			s = s.substring(s.lastIndexOf(":", i) + 1, i);

			return HGNC.getSymbol(s);
		}
		return null;
	}
}
