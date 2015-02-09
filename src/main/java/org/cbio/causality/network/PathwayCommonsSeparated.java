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
public class PathwayCommonsSeparated
{
	private static final String url = "http://cbio.mskcc.org/~ozgun/PCWithLoc.sif.zip";
	private static final String tempFile = "PCWithLoc.sif";
	private static final String dir = "PC-res/";

	public static Graph getGraph(Set<String> sources, SIFType... types)
	{
		if (fileExists(types))
		{
			if (types.length == 1) return getSingleGraph(types[0], sources);
			else if (types.length > 1)
			{
				GraphList graph = new GraphList("Pathway Commons");

				for (SIFType type : types)
				{
					graph.addGraph(getSingleGraph(type, sources));
				}

				return graph;
			}
		}
		return null;
	}

	public static Graph getSingleGraph(SIFType type, Set<String> sources)
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

			boolean valid = false;
			for (String res : token[2].split(";"))
			{
				if (sources.contains(res))
				{
					valid = true;
					break;
				}
			}

			if (!valid) continue;

			graph.putRelation(token[0], token[1], type.isDirected());
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

			Set<String> res = new HashSet<String>();

			while (sc.hasNextLine())
			{
				String line = sc.nextLine();
				String[] token = line.split("\t");

				if (token.length > 2)
				{
					if (!writers.containsKey(token[1])) writers.put(token[1],
						new BufferedWriter(new FileWriter(dir + token[1] + ".txt")));

					writers.get(token[1]).write(token[0] + "\t" + token[2]);

					if (token.length > 5)
					{
						writers.get(token[1]).write("\t" + token[5]);
						Collections.addAll(res, token[5].split(";"));
					}

					writers.get(token[1]).write("\n");
				}
			}

			for (Writer writer : writers.values())
			{
				writer.close();
			}

			new File(tempFile).delete();

			for (String re : res)
			{
				System.out.println("resource = " + re);
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return false;
	}

	public static void main(String[] args)
	{
		printNetworkSizes();
	}

	private static void printNetworkSizes()
	{
		Graph graph = getGraph(new HashSet<String>(Arrays.asList("reactome")), SIFEnum.values());
		graph.printStats();

		graph = getGraph(new HashSet<String>(Arrays.asList("humancyc")), SIFEnum.values());
		graph.printStats();
	}
}
