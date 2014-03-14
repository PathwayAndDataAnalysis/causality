package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
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

	public static Graph getGraph(SIFType... types)
	{
		String edgeType = types[0].getTag();

		for (int i = 1; i < types.length; i++)
		{
			edgeType += "," + types[i].getTag();
		}

		if (fileExists(types))
		{
			Graph graph = new Graph("Signed PC", edgeType);

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

					graph.putRelation(token[0], token[1], type.isDirected());
				}
			}

			return graph;
		}
		return null;
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

					writers.get(token[1]).write(token[0] + "\t" + token[2] + "\n");
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
