package org.cbio.causality.data;

import org.cbio.causality.data.portal.CBioPortalManager;
import org.cbio.causality.util.Waiter;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * Beware: This class is not tread-safe.
 * @author Ozgun Babur
 */
public class GeneCards
{
	private static final Set<String> CANCER_KW = new HashSet<String>(Arrays.asList(
		"cancer", "oma", "oplasm", "tumor"));

	private static final String NOT_ASSOCIATED_FILE = "not-associated.txt";

	public static String getCacheDir()
	{
		return CBioPortalManager.getCacheDir() + File.separator + "GeneCards" + File.separator;
	}

	public static Set<String> getRelatedCancers(String gene)
	{
		Set<String> cancers = readCancersFromCache(gene);
		if (cancers != null && cancers.isEmpty()) return cancers;

		if (cancers == null)
		{
			try
			{
				cancers = spiderCancers(gene);
			}
			catch (RuntimeException e)
			{
				e.printStackTrace();
			}
		}

		if (cancers != null)
		{
			if (cancers.isEmpty()) writeNotAssociated(gene);
			else cacheCancers(gene, cancers);
		}

		return cancers;
	}



	private static Set<String> readCancersFromCache(String gene)
	{
		Set<String> na = getNotAssociated();
		if (na.contains(gene)) return Collections.emptySet();
		else
		{
			String filename = getCacheDir() + gene + ".txt";
			if (new File(filename).exists())
			{
				try
				{
					Set<String> set = new HashSet<String>();
					BufferedReader reader = new BufferedReader(new FileReader(filename));

					for (String line = reader.readLine(); line != null; line = reader.readLine())
					{
						if (!line.isEmpty()) set.add(line);
					}

					reader.close();
					return set;
				}
				catch (IOException e)
				{
					e.printStackTrace();
					return null;
				}
			}
			else return null;
		}
	}

	private static void cacheCancers(String gene, Set<String> cancers)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(
				new FileWriter(getCacheDir() + gene + ".txt"));

			for (String cancer : cancers)
			{
				writer.write(cancer + "\n");
			}

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private static Set<String> getNotAssociated()
	{
		Set<String> genes = new HashSet<String>();

		try
		{
			String fileName = getCacheDir() + NOT_ASSOCIATED_FILE;

			if (new File(fileName).exists())
			{

				BufferedReader reader = new BufferedReader(
					new FileReader(fileName));

				for (String line = reader.readLine(); line != null; line = reader.readLine())
				{
					if (!line.isEmpty()) genes.add(line);
				}

				reader.close();
			}
			else
			{
				BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));

				writer.close();
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return genes;
	}

	private static void writeNotAssociated(String gene)
	{
		Set<String> set = getNotAssociated();
		set.add(gene);
		try
		{
			BufferedWriter writer = new BufferedWriter(
				new FileWriter(getCacheDir() + NOT_ASSOCIATED_FILE));

			for (String s : set)
			{
				writer.write(s + "\n");
			}

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private static Set<String> spiderCancers(String gene)
	{
		Waiter.pause(5000);
		Set<String> cancers = new HashSet<String>();
		try
		{
			URL url = new URL("http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene);
			HttpURLConnection con = (HttpURLConnection) url.openConnection();

			con.addRequestProperty("User-Agent", "Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10.4; en-US; rv:1.9.2.2) Gecko/20100316 Firefox/3.6.2");

			BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()));

			StringBuilder content = null;
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				if (content == null && line.equals("<br ><h2 class=\"navbar\">Disorders "))
				{
					content = new StringBuilder();
				}
				else if (content != null)
				{
					if (line.equals("\t<TR>"))
					{
						cancers = parseCancers(content.toString());
					}
					else content.append(line);
				}
			}

			reader.close();
		}
		catch (IOException e)
		{
			throw new RuntimeException("Error during spidering.", e);
		}
		return cancers;
	}

	private static Set<String> parseCancers(String line)
	{
		Set<String> words = new HashSet<String>();

		int i = 0;

		while (i >= 0 && i < line.length() - 1)
		{
			i = line.indexOf(">", i + 1);
			if (i >= 0)
			{
				int j = line.indexOf("<", i + 1);

				if (j > 0)
				{
					String word = line.substring(i + 1, j);
					word = word.toLowerCase();

					for (String kw : CANCER_KW)
					{
						if (word.endsWith(kw))
						{
							words.add(word);
							break;
						}
					}
				}
			}
		}

		return words;
	}

	static
	{
		File dir = new File(getCacheDir());
		if (!dir.exists()) dir.mkdirs();
	}

	public static void main(String[] args)
	{
		String[] genes = "DTNA DYTN UTRN DTNB DRP2".split(" ");


		for (String gene : genes)
		{
			if (gene.isEmpty()) continue;

			System.out.print("\n" + gene + ":");
			for (String disease : getRelatedCancers(gene))
			{
				System.out.print(" " + disease + ",");
			}
		}
	}
}
