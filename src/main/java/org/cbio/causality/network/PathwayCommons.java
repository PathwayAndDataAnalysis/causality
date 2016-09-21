package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.GraphList;
import org.cbio.causality.util.Download;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PathwayCommons
{
	private static final String url = "http://www.pathwaycommons.org/archives/PC2/current/PathwayCommons.8.All.EXTENDED_BINARY_SIF.hgnc.txt.gz";
	private static final String tempFile = "PC.sif";
	private static final String dir = "PC/";

	public static Graph getGraph(SIFType... types)
	{
		if (fileExists(types))
		{
			if (types.length == 1) return getSingleGraph(types[0]);
			else if (types.length > 1)
			{
				GraphList graph = new GraphList("Pathway Commons");

				for (SIFType type : types)
				{
					graph.addGraph(getSingleGraph(type));
				}

				return graph;
			}
		}
		return null;
	}

	public static Graph getSingleGraph(SIFType type)
	{
		String edgeType = type.getTag();

		Graph graph = new Graph("Pathway Commons", edgeType);

		Scanner sc;
		try
		{
			sc = new Scanner(new FileInputStream(dir + type.getTag() + ".txt"));
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

			if (token.length > 2)
			{
				graph.putRelation(token[0], token[1], token[2], type.isDirected());
			}
			else
			{
				graph.putRelation(token[0], token[1], type.isDirected());
			}
		}

		return graph;
	}

	private static boolean fileExists(SIFType... types)
	{
		if (!fileExistsJustCheck(types))
		{
			if (!(new File(tempFile).exists())) Download.downloadAndUncompress(url, tempFile);
			extractData();
		}

		return fileExistsJustCheck(types);
	}

	private static boolean fileExistsJustCheck(SIFType[] types)
	{
		boolean present = true;
		for (SIFType type : types)
		{
			if (!new File(dir + type.getTag() + ".txt").exists())
			{
				present = false;
				break;
			}
		}
		return present;
	}

	private static boolean extractData()
	{
		try
		{
			Scanner sc = new Scanner(new File(tempFile));

			Map<String, Writer> writers = new HashMap<String, Writer>();

			new File(dir).mkdirs();

			while (sc.hasNextLine())
			{
				String line = sc.nextLine();
				if (line.isEmpty()) break;

				String[] token = line.split("\t");

				if (token.length > 2)
				{
					if (!writers.containsKey(token[1])) writers.put(token[1],
						new BufferedWriter(new FileWriter(dir + token[1] + ".txt")));

					writers.get(token[1]).write(token[0] + "\t" + token[2]);

					if (token.length > 6)
					{
						writers.get(token[1]).write("\t" + token[6]);
					}

					writers.get(token[1]).write("\n");
				}
			}

			for (Writer writer : writers.values())
			{
				writer.close();
			}

			new File(tempFile).delete();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return false;
	}

	public static void main(String[] args)
	{
//		printDataOverlaps();
		printNetworkSizes();
//		printMostConnected(getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF), 500);
	}

	private static void printNetworkSizes()
	{
		Graph graph = getGraph(SIFEnum.values());
		graph.printStats();
	}

	private static void printDataOverlaps()
	{
		Graph graph = getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		graph.printVennIntersections(true, SPIKE.getGraphPostTl(), SignaLink.getGraphPostTl(), ReactomeFI.getGraphPostTl());

		System.out.println();
		graph = getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		graph.printVennIntersections(true, SPIKE.getGraphTR(), SignaLink.getGraphPostTl(), ReactomeFI.getGraphTR(), MSigDBTFT.getGraph());

		System.out.println();
		Graph icw = PathwayCommons.getGraph(SIFEnum.IN_COMPLEX_WITH);
		icw.printVennIntersections(false, HPRD.getGraph(), IntAct.getGraph(), ReactomeFI.getGraphPPI());
	}

	private static void printMostConnected(Graph graph, int limit)
	{
		List<String> genes = new ArrayList<String>(graph.getSymbols());
		final Map<String, Integer> degree = new HashMap<String, Integer>();
		for (String gene : genes)
		{
			degree.put(gene, graph.getDegree(gene));
		}
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return degree.get(o2).compareTo(degree.get(o1));
			}
		});

		int i = 0;
		for (String gene : genes)
		{
			i++;
			System.out.println(gene + "\t" + degree.get(gene));
			if (i == limit) break;
		}
	}
}
