package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SignaLink
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
		graphPostTl = new Graph("SignaLink post-translational mod", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graphTR = new Graph("SignaLink transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream("signalink.txt"));

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

	//--- Section: Prepare resources --------------------------------------------------------------|

	private static void prepare() throws IOException
	{
		Scanner sc = new Scanner(SignaLink.class.getResourceAsStream("signalink-raw.txt"));
		BufferedWriter writer = new BufferedWriter(new FileWriter("signalink.txt"));

		Set<String> existing = new HashSet<String>();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split(";");

			String source = HGNC.getSymbol(token[0]);
			String target = HGNC.getSymbol(token[6]);

			if (source == null || target == null) continue;

			String type = token[13].toLowerCase();
			if (type.contains("undirected")) continue;

			boolean tr = type.contains("transcriptional");

			String rel = source + "\t" + target + "\t" + (tr ? "T" : "");

			if (existing.contains(rel)) continue;
			else existing.add(rel);

			writer.write(rel + (sc.hasNextLine() ? "\n" : ""));
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		System.out.println(getGraphPostTl().getDownstream("PPM1E"));
	}
}
