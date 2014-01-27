package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

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
		graph = new Graph("SPIKE", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream("SPIKE-ParsedFromXML.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[1])) continue;

			token[0] = HGNC.getSymbol(token[0]);
			token[1] = HGNC.getSymbol(token[1]);

			if (token[0] == null || token[1] == null) continue;

			graph.putRelation(token[0], token[1], true);
		}
	}
}
