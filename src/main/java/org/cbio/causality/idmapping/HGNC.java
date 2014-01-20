package org.cbio.causality.idmapping;

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
	private static Map<String, String> id2sym;
	private static Map<String, String> old2new;
	private static Map<String, Set<String>> families = new HashMap<String, Set<String>>();

	/**
	 * Gets the latest approved official symbol related to the given ID or symbol. If the parameter
	 * is ID, then it should start with "HGNC:".
	 * @param symbolOrID HGNC ID, symbol, or a previous symbol
	 * @return latest symbol
	 */
	public static String getSymbol(String symbolOrID)
	{
		if (id2sym.containsKey(symbolOrID)) return id2sym.get(symbolOrID);
		else if (sym2id.containsKey(symbolOrID)) return symbolOrID;
		symbolOrID = symbolOrID.toUpperCase();
		if (old2new.containsKey(symbolOrID)) return old2new.get(symbolOrID);
		return null;
	}

	public static Set<String> getFamily(String name)
	{
		if (!families.containsKey(name)) return Collections.emptySet();
		return families.get(name);
	}

	static
	{
		try
		{
			sym2id = new HashMap<String, String>();
			id2sym = new HashMap<String, String>();
			old2new = new HashMap<String, String>();
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
						old2new.put(old, sym);
					}
				}
				if (token.length > 3)
				{
					for (String synonym : token[3].split(","))
					{
						synonym = synonym.trim().toUpperCase();
						old2new.put(synonym, sym);
					}
				}
				if (token.length > 4)
				{
					if (!families.containsKey(token[4]))
						families.put(token[4], new HashSet<String>());

					families.get(token[4]).add(sym);
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

	public static void main(String[] args)
	{
		System.out.println(HGNC.getSymbol("P53"));
	}
}
