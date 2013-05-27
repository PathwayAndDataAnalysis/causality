package org.cbio.causality.hprd;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ShuffledHPRD implements InteractionProvider
{
	private Map<String, String> redirectFwd;
	private Map<String, String> redirectBkw;
	private Map<String, Set<String>> map;

	public ShuffledHPRD()
	{
		initShuffledPreserveDegrees();
	}

	private void initShuffledPlain()
	{
		List<String> list1 = new ArrayList<String>(HPRD.getAllSymbols());
		List<String> list2 = new ArrayList<String>(list1);
		Collections.shuffle(list2);

		redirectFwd = new HashMap<String, String>(list1.size());
		redirectBkw = new HashMap<String, String>(list1.size());

		for (int i = 0; i < list1.size(); i++)
		{
			redirectFwd.put(list1.get(i), list2.get(i));
			redirectBkw.put(list2.get(i), list1.get(i));
		}

		map = new HashMap<String, Set<String>>(list1.size());
	}

	private void initShuffledPreserveDegrees()
	{
		List<String> genes = new ArrayList<String>(HPRD.getAllSymbols());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return ((Integer) HPRD.map.get(o2).size()).compareTo(HPRD.map.get(o1).size());
			}
		});

		int size = genes.size();
		int bins = 200;
		List<String>[] list = new ArrayList[bins];

		for (int i = 0; i < bins; i++)
		{
			list[i] = new ArrayList<String>(genes.subList(
				(int) Math.floor(i * size / (double) bins),
				(int) Math.floor((i + 1) * size / (double) bins)));

			Collections.shuffle(list[i]);
		}

		List<String> shuffled = new ArrayList<String>();
		for (List<String> aList : list)	shuffled.addAll(aList);

		redirectFwd = new HashMap<String, String>(genes.size());
		redirectBkw = new HashMap<String, String>(genes.size());

		for (int i = 0; i < genes.size(); i++)
		{
			redirectFwd.put(genes.get(i), shuffled.get(i));
			redirectBkw.put(shuffled.get(i), genes.get(i));
		}

		map = new HashMap<String, Set<String>>(genes.size());
	}

	public Set<String> getInteractions(String symbol)
	{
		if (map.containsKey(symbol)) return new HashSet<String>(map.get(symbol));

		String s = redirectBkw.get(symbol);

		Set<String> inter = HPRD.getInteractors(s);

		Set<String> redirInter = new HashSet<String>(inter.size());

		for (String s1 : inter)
		{
			redirInter.add(redirectFwd.get(s1));
		}
		map.put(symbol, redirInter);

		return redirInter;
	}
}
