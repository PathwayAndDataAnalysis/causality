package org.cbio.causality.idmapping;

import org.cbio.causality.analysis.Graph;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class HGNC
{
	private static Map<String, String> sym2id;
	private static Map<String, String> sym2chr;
	private static Map<String, String> id2sym;
	private static Map<String, String> old2new;
	private static Map<String, String> uniprot2sym;
	private static Map<String, Set<String>> families = new HashMap<String, Set<String>>();


	/**
	 * Gets the latest approved official symbol related to the given ID or symbol. If the parameter
	 * is ID, then it should start with "HGNC:".
	 * @param symbolOrID HGNC ID, symbol, or a previous symbol
	 * @return latest symbol
	 */
	public static String getSymbol(String symbolOrID)
	{
		if (symbolOrID == null) return null;
		if (id2sym.containsKey(symbolOrID)) return id2sym.get(symbolOrID);
		else if (sym2id.containsKey(symbolOrID)) return symbolOrID;
		symbolOrID = symbolOrID.toUpperCase();
		if (old2new.containsKey(symbolOrID)) return old2new.get(symbolOrID);
		if (uniprot2sym.containsKey(symbolOrID)) return uniprot2sym.get(symbolOrID);
		return null;
	}

	public static Set<String> getFamily(String name)
	{
		if (!families.containsKey(name)) return Collections.emptySet();
		return families.get(name);
	}

	public static Set<String> getAllSymbols()
	{
		return sym2id.keySet();
	}

	public static Set<String> getSymbolsOfChromosome(String chr)
	{
		Set<String> set = new HashSet<String>();
		for (String sym : getAllSymbols())
		{
			if (sym2chr.containsKey(sym))
			{
				String c = sym2chr.get(sym);
				String no = c.split("p")[0].split("q")[0];

				String desNo = chr.split("p")[0].split("q")[0];

				if (no.equals(desNo) && c.contains(chr))
				{
					set.add(sym);
				}
			}
		}
		return set;
	}

	public static String getChromosomeLoc(String symbol)
	{
		return sym2chr.get(symbol);
	}

	static
	{
		try
		{
			sym2id = new HashMap<String, String>();
			sym2chr = new HashMap<String, String>();
			id2sym = new HashMap<String, String>();
			old2new = new HashMap<String, String>();
			uniprot2sym = new HashMap<String, String>();
			families = new HashMap<String, Set<String>>();

			BufferedReader reader = new BufferedReader(new InputStreamReader(
				HGNC.class.getResourceAsStream("hgnc.txt")));

			reader.readLine(); //skip header
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				String[] token = line.split("\t");
				String sym = token[1];
				String id = token[0];
				sym2id.put(sym, id);
				id2sym.put(id, sym);

				if (token.length > 2)
				{
					for (String old : token[2].split(","))
					{
						old = old.trim().toUpperCase();
						if (old.isEmpty()) continue;
						old2new.put(old, sym);
					}
				}
				if (token.length > 3)
				{
					for (String synonym : token[3].split(","))
					{
						synonym = synonym.trim().toUpperCase();
						if (synonym.isEmpty()) continue;
						old2new.put(synonym, sym);
					}
				}
				if (token.length > 4)
				{
					sym2chr.put(sym, token[4]);
				}
				if (token.length > 5)
				{
					if (!families.containsKey(token[5]))
						families.put(token[5], new HashSet<String>());

					families.get(token[5]).add(sym);
				}
				if (token.length > 6 && !token[6].isEmpty())
				{
					uniprot2sym.put(token[6], sym);
				}
			}
			reader.close();
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public static Graph getCompleteClique(boolean directed)
	{
		Graph graph = new Graph("HGNC complete clique", "clique-edge");
		for (String s1 : sym2id.keySet())
		{
			for (String s2 : sym2id.keySet())
			{
				if (s1.equals(s2)) continue;
				if (!directed && s2.compareTo(s1) < 0) continue;
				graph.putRelation(s1, s2, directed);
				graph.putRelation(s1, s2, directed);
			}
		}
		return graph;
	}

	public static void main(String[] args)
	{
//		System.out.println(getSymbol("CDC42"));
//		System.out.println(getChromosomeLoc("ATAD2"));
		Set<String> set = getSymbolsOfChromosome("8q24");
//		for (String sym : set)
//		{
//			System.out.println(sym + "\t" + getChromosomeLoc(sym));
//		}

		System.out.println(getChromosomeLoc("MYC"));
	}
}
