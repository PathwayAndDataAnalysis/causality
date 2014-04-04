package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.PhosphoGraph;
import org.cbio.causality.signednetwork.SignedType;
import org.cbio.causality.util.Download;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class SignedPC extends PathwayCommons
{
	private static final String url = "http://resources.chibe.googlecode.com/hg/SignedPC.sif.gz";
	private static final String tempFile = "SignedPC.sif";
	private static final String dir = "SignedPC/";

	public static Map<SignedType, Graph> getAllGraphs()
	{
		Map<SignedType, Graph> map = new HashMap<SignedType, Graph>();
		for (SignedType type : SignedType.values())
		{
			map.put(type, getGraph(type));
		}
		return map;
	}

	public static Graph getGraph(SIFType... types)
	{
		String edgeType = types[0].getTag();

		for (int i = 1; i < types.length; i++)
		{
			edgeType += "," + types[i].getTag();
		}

		boolean phos = false;
		for (SIFType type : types)
		{
			if (type == SignedType.PHOSPHORYLATES || type == SignedType.DEPHOSPHORYLATES)
			{
				phos = true;
				break;
			}
		}
		
		if (fileExists(types))
		{
			Graph graph = phos ? new PhosphoGraph("Signed PC", edgeType) :
				new Graph("Signed PC", edgeType);

			for (SIFType type : types)
			{
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
						if (phos && token.length > 3)
						{
							((PhosphoGraph) graph).putRelation(token[0], token[1], token[2], type.isDirected(), token[3]);
						}
						else
						{
							graph.putRelation(token[0], token[1], token[2], type.isDirected());
						}
					}
					else
					{
						graph.putRelation(token[0], token[1], type.isDirected());
					}
				}
			}

			return graph;
		}
		return null;
	}

	public static Graph getSingleGraph(SIFType type)
	{
		return getGraph(type);
	}

	private static boolean fileExists(SIFType[] types)
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
				String[] token = line.split("\t");

				if (token.length > 2)
				{
					if (!writers.containsKey(token[1])) writers.put(token[1],
						new BufferedWriter(new FileWriter(dir + token[1] + ".txt")));

					writers.get(token[1]).write(token[0] + "\t" + token[2]);

					if (token.length > 3)
					{
						writers.get(token[1]).write("\t" + token[3]);
					}
					if (token.length > 4)
					{
						writers.get(token[1]).write("\t" + token[4]);
					}

					writers.get(token[1]).write("\n");
				}
			}

			for (Writer writer : writers.values())
			{
				writer.close();
			}

//			new File(tempFile).delete();
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
		for (SignedType type : SignedType.values())
		{
			getGraph(type).printStats();
		}
	}
}
