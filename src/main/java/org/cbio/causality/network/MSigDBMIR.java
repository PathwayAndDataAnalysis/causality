package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.util.Scanner;
import java.util.Set;

/**
 * MicroRNA targets.
 * @author Ozgun Babur
 */
public class MSigDBMIR
{
	private static Graph graph;

	public static Graph getGraph()
	{
		return graph;
	}

	static
	{
		graph = new Graph("MIR", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		Scanner sc = new Scanner(MSigDBMIR.class.getResourceAsStream("msigdb-mir.gmt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			for (String mir : token[0].split(","))
			{
				if (mir.startsWith("MIR-"))
				{
					for (int i = 2; i < token.length; i++)
					{
						token[i] = HGNC.getSymbol(token[i]);
						if (token[i] == null) continue;

						graph.putRelation(mir, token[i], true);
					}
				}
			}
		}
	}

	public static void main(String[] args)
	{
		Graph graph = getGraph();
		System.out.println(graph.getDownstream("MIR-17-3P"));
	}
}
