package org.cbio.causality.network;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.util.Scanner;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class MSigDBTFT
{
	private static Graph graph;

	public static Graph getGraph()
	{
		return graph;
	}

	static
	{
		graph = new Graph();

		Scanner sc = new Scanner(MSigDBTFT.class.getResourceAsStream("msigdb-tft.gmt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			String s = token[0];

			int dInd = s.indexOf("$");
			if (dInd < 0) continue;
			int uInd = s.indexOf("_", dInd);
			if (uInd < 0) continue;

			s = s.substring(dInd + 1, uInd);

			String tf = HGNC.getSymbol(s);

			if (tf != null)
			{
				for (int i = 2; i < token.length; i++)
				{
					addRelation(token[i], tf);
				}
			}
			else
			{
				Set<String> family = HGNC.getFamily(s);

				for (String mem : family)
				{
					for (int i = 2; i < token.length; i++)
					{
						addRelation(token[i], mem);
					}
				}
			}
		}
	}

	private static void addRelation(String target, String tf)
	{
		target = HGNC.getSymbol(target);
		if (target == null) return;
		graph.putRelation(tf, target, true);
	}

	public static void main(String[] args)
	{
		Graph graph = getGraph();
		System.out.println(graph.getDownstream("MEF2C"));
	}
}
