package org.cbio.causality.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.Download;
import org.cbio.causality.util.FileUtil;

import java.io.*;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class ReactomeFI
{
	public static final String resourceFile = "ReactomeFI.txt";
	public static final String tempFile = "ReactomeFI.zip";
	public static final String url = "http://reactomews.oicr.on.ca:8080/caBigR3WebApp2013/FIsInGene_121013_with_annotations.txt.zip";

	private static Graph graphPostTl;
	private static Graph graphTR;
	private static Graph graphPPI;

	public static Graph getGraphPostTl()
	{
		return graphPostTl;
	}

	public static Graph getGraphTR()
	{
		return graphTR;
	}

	public static Graph getGraphPPI()
	{
		return graphPPI;
	}

	static
	{
		File file = new File(resourceFile);
		if (!file.exists())
		{
			if (Download.downloadAsIs(url, tempFile))
			{
				FileUtil.extractEntryContainingNameInZipFile(tempFile, "FI", "MACOSX", resourceFile);
			}
		}

		try
		{
			Scanner sc = new Scanner(file);

			graphPostTl = new Graph("Reactome FI PPrel", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
			graphTR = new Graph("Reactome FI GErel", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());
			graphPPI = new Graph("Reactome FI PPI", SIFEnum.INTERACTS_WITH.getTag());

			while (sc.hasNextLine())
			{
				String line = sc.nextLine();

				String[] token = line.split("\t");

				if (token[0].equals(token[1])) continue;

				token[0] = HGNC.getSymbol(token[0]);
				token[1] = HGNC.getSymbol(token[1]);

				if (token[0] == null || token[1] == null) continue;

				String arrow = token[3];
				String annot = token[2];

				if (arrow.contains("->") || arrow.contains("-|"))
				{
					if (annot.contains("expression regulates") ||
						annot.contains("GErel: expression"))
					{
						graphTR.putRelation(token[0], token[1], true);
					}
					else
					{
						graphPostTl.putRelation(token[0], token[1], true);
					}
				}
				if (arrow.contains("<-") || arrow.contains("|-"))
				{
					if (annot.contains("expression regulated by") ||
						annot.contains("GErel: expression by"))
					{
						graphTR.putRelation(token[1], token[0], true);
					}
					else
					{
						graphPostTl.putRelation(token[1], token[0], true);
					}
				}
				if (arrow.equals("-"))
				{
					graphPPI.putRelation(token[0], token[1], false);
				}
			}
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
	}

	//--- Section: Prepare resources --------------------------------------------------------------|

	public static void main(String[] args) throws IOException
	{
		ReactomeFI.getGraphPostTl().printStats();
		System.out.println("--------");
		ReactomeFI.getGraphTR().printStats();
	}
}
