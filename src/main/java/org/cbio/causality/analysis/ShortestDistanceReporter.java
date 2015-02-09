package org.cbio.causality.analysis;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ShortestDistanceReporter
{
	public static Map<String, Map<String, Integer>> getShortestDistances(Graph graph, int limit)
	{
		Map<String, Map<String, Integer>> dist = new HashMap<String, Map<String, Integer>>();

		for (String source : graph.getSymbols())
		{
			for (int i = 1; i <= limit; i++)
			{
				Set<String> dw = graph.getDownstream(Collections.singleton(source), i);

				for (String target : dw)
				{
					if (!dist.containsKey(source)) dist.put(source, new HashMap<String, Integer>());
					if (!dist.get(source).containsKey(target)) dist.get(source).put(target, i);
				}
			}
		}
		return dist;
	}

	public static void writeShortestDistances(Graph graph, int limit, String filename) throws IOException
	{
		Map<String, Map<String, Integer>> dist = getShortestDistances(graph, limit);
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (String source : dist.keySet())
		{
			for (String target : dist.get(source).keySet())
			{
				writer.write(source + "\t" + target + "\t" + dist.get(source).get(target) + "\n");
			}
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		Graph graph = new Graph("PERA", "is-upstream-of");
		Scanner sc = new Scanner(new File("/home/ozgun/Downloads/prior_network.tsv"));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");
			graph.putRelation(token[0], token[2], true);
		}

		graph.printStats();

		writeShortestDistances(graph, 10, "/home/ozgun/Downloads/distances.txt");
	}
}
