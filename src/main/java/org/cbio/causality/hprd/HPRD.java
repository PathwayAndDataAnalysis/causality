package org.cbio.causality.hprd;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class HPRD implements InteractionProvider
{
	protected static Map<String, Set<String>> map;

	public static Set<String> getInteractors(String symbol)
	{
		if (map.containsKey(symbol)) return new HashSet<String>(map.get(symbol));
		else return Collections.emptySet();
	}

	public static Set<String> getAllSymbols()
	{
		return map.keySet();
	}

	public static void cropTo(Set<String> genes)
	{
		for (String s : new HashSet<String>(map.keySet()))
		{
			if (!genes.contains(s)) map.remove(s);
		}
		for (String s : map.keySet())
		{
			map.get(s).retainAll(genes);
		}
		for (String s : new HashSet<String>(map.keySet()))
		{
			if (map.get(s).isEmpty()) map.remove(s);
		}
	}

	public static int getDegree(String symbol)
	{
		if (map.containsKey(symbol)) return map.get(symbol).size();
		return 0;
	}

	static
	{
		map = new HashMap<String, Set<String>>();

		Scanner sc = new Scanner(HPRD.class.getResourceAsStream(
			"BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			String[] token = line.split("\t");

			if (token[0].equals(token[3])) continue;

			if (!map.containsKey(token[0])) map.put(token[0], new HashSet<String>());
			if (!map.containsKey(token[3])) map.put(token[3], new HashSet<String>());

			map.get(token[0]).add(token[3]);
			map.get(token[3]).add(token[0]);
		}
	}

	@Override
	public Set<String> getInteractions(String symbol)
	{
		return HPRD.getInteractors(symbol);
	}

	public static int getEdgeSize()
	{
		int size = 0;
		for (Set<String> set : map.values())
		{
			size += set.size();
		}
		assert size % 2 == 0;
		return size / 2;
	}

	public static void main(String[] args)
	{
		Set<String> set = HPRD.getInteractors("BRD4");
		for (String s : set)
		{
			System.out.println("s = " + s);
		}
	}
}
