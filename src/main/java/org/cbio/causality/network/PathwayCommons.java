package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.util.Download;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class PathwayCommons
{
	private static final String url = "http://resources.chibe.googlecode.com/hg/PC.sif.gz";
	private static final String filename = "PC.sif";

	public static Graph getGraph(SIFType... types)
	{
		String edgeType = types[0].getTag();

		for (int i = 1; i < types.length; i++)
		{
			edgeType += "," + types[i].getTag();
		}

		if (fileExists())
		{
			Graph graph = new Graph("Pathway Commons", edgeType);
			Set<SIFType> typeSet = new HashSet<SIFType>(Arrays.asList(types));

			Scanner sc;
			try
			{
				sc = new Scanner(new FileInputStream(filename));
			}
			catch (FileNotFoundException e)
			{
				e.printStackTrace();
				return null;
			}
			while (sc.hasNextLine())
			{
				String line = sc.nextLine();

				String[] token = line.split("\t");

				SIFType type = SIFEnum.typeOf(token[1]);

				if (typeSet.contains(type))
				{
					graph.putRelation(token[0], token[2], type.isDirected());
				}
			}

			return graph;
		}
		return null;
	}

	private static boolean fileExists()
	{
		return new File(filename).exists() || Download.downloadAndUncompress(url, filename);
	}

	public static void main(String[] args)
	{
		Graph graph = getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF);
		System.out.println("graph.getNeighbors(\"TP53\").size() = " + graph.getNeighbors("TP53").size());
	}
}
