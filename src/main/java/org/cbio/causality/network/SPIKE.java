package org.cbio.causality.network;

import org.cbio.causality.analysis.Graph;

import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class SPIKE
{
	private static Graph graph;

	public static Graph getGraph()
	{
		return graph;
	}

	static
	{
		graph = new Graph();

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream("SPIKE-ParsedFromXML.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[1])) continue;

			graph.putRelation(token[0], token[1], true);
		}
	}
}
