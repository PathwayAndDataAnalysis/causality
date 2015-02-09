package org.cbio.causality.signednetwork;

import org.cbio.causality.idmapping.CancerGeneCensus;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * I am implementing this class because Chris wants to know the related paper references for
 * causality relations.
 *
 * @author Ozgun Babur
 */
public class PubmedExtractor
{
	public static void main(String[] args) throws IOException
	{
		printReferences("/home/ozgun/Downloads/liposarcoma/20150108_LiposarcomaXenografts_v2.sif",
			"SignedPC.sif");
	}

	public static void printReferences(String subsetFile, String signedPCFileWithRefs) throws IOException
	{
		Set<String> subset = readSubset(subsetFile);
		System.out.println("subset.size() = " + subset.size());
		final Map<String, Set<String>> map = getPubmed2Rels(signedPCFileWithRefs, subset);

		List<String> pmids = new ArrayList<String>(map.keySet());
		Collections.sort(pmids, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Integer(map.get(o2).size()).compareTo(map.get(o1).size());
			}
		});

		Set<String> covered = new HashSet<String>();
		for (Set<String> rels : map.values())
		{
			covered.addAll(rels);
		}

		System.out.println("Not covered:");
		for (String rel : subset)
		{
			if (!covered.contains(rel))
			{
				System.out.println(rel.replaceAll("\t", " "));
			}
		}
		System.out.println();

//		for (String pmid : pmids)
//		{
//			System.out.println(getArticleDetails(pmid));
//			System.out.println("Covered interactions:");
//			for (String rel : map.get(pmid))
//			{
//				System.out.println(rel.replaceAll("\t", " "));
//			}
//			System.out.println("\n");
//		}
	}

	private static Set<String> readSubset(String file) throws FileNotFoundException
	{
		Set<String> set = new HashSet<String>();
		Scanner sc = new Scanner(new File(file));
		Set<String> genes = new HashSet<String>();
		
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			if (token.length > 2)
			{
				set.add(token[0] + "\t" + token[1] + "\t" + token[2]);
				genes.add(token[0]);
				genes.add(token[2]);
			}
			else if (token.length == 1)
			{
//				genes.add(token[0]);
			}
		}

		Set<String> census = CancerGeneCensus.getAllSymbols();
		genes.retainAll(census);
		System.out.println(new ArrayList<String>(genes));

		return set;
	}
	
	private static Map<String, Set<String>> getPubmed2Rels(String file, Set<String> subset) throws FileNotFoundException
	{
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();
		Scanner sc = new Scanner(new File(file));

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			if (token.length > 5)
			{
				String id = token[0] + "\t" + token[1] + "\t" + token[2];
				
				if (subset.contains(id))
				{
					for (String pmid : token[5].split(" "))
					{
						if (!map.containsKey(pmid)) map.put(pmid, new HashSet<String>());
						map.get(pmid).add(id);
					}
				}
			}
		}
		
		return map;
	}

	private static String getArticleDetails(String pmid) throws IOException
	{
		URL url = new URL("http://www.ncbi.nlm.nih.gov/pubmed/" + pmid);
		URLConnection con = url.openConnection();
		Scanner sc = new Scanner(con.getInputStream());

		String title = null;
		String trunc = null;
		String auth = null;
		String details = null;

		String authprefix = "<meta name=\"author\" content=\"";
		String detprefix = "<meta name=\"description\" content=\"";

		String[] unwantedEnds = new String[]{" Research Support", " Comparative Study\"",
			" Comparative Study;", " Evaluation Studies;", " Epub ", " doi: ", };

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			if (line.startsWith("    <title>"))
			{
				title = line.substring(line.indexOf(">") + 1,
					line.length() - "  - PubMed - NCBI".length());

				if (title.endsWith("...")) trunc = title.substring(0, title.length() - 3);
			}
			else if (trunc != null && line.contains(trunc))
			{
				title = line.substring(line.indexOf(trunc), line.indexOf("<", line.indexOf(trunc)));
				trunc = null;
			}
			else if (line.contains(authprefix))
			{
				auth = line.substring(line.indexOf(authprefix) + authprefix.length(), line.indexOf("\"", line.indexOf(authprefix) + authprefix.length()));
			}

			if (line.contains(detprefix))
			{
				details = line.substring(line.indexOf(detprefix) + detprefix.length(), line.indexOf("\"", line.indexOf(detprefix) + detprefix.length()));
				for (String s : unwantedEnds)
				{
					if (details.contains(s)) details = details.substring(0, details.indexOf(s));
				}
			}


			if (title != null && auth != null && details != null && trunc == null) break;
		}

		String s = title + "\n" + auth + "\n" + details + "\n" + url;
		return s;
	}
}
