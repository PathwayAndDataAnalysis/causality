package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.io.IOException;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class SPIKE
{
	private static Graph graphPostTl;
	private static Graph graphTR;

	public static Graph getGraphPostTl()
	{
		return graphPostTl;
	}

	public static Graph getGraphTR()
	{
		return graphTR;
	}

	static
	{
		graphPostTl = new Graph("SPIKE Post-translational mod", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graphTR = new Graph("SPIKE transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream("SPIKE-ParsedFromXML.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[1])) continue;

			token[0] = HGNC.getSymbol(token[0]);
			token[1] = HGNC.getSymbol(token[1]);

			if (token[0] == null || token[1] == null) continue;

			boolean tr = token.length > 2 && token[2].equals("T");

			if (tr) graphTR.putRelation(token[0], token[1], true);
			else graphPostTl.putRelation(token[0], token[1], true);
		}
	}

	public static void main(String[] args) throws IOException
	{
		System.out.println(getGraphPostTl().getDownstream("PPM1E"));
	}
}
