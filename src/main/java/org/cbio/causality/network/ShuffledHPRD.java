package org.cbio.causality.network;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
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

		initShuffledMapping(list1, list2);
	}

	private void initShuffledPreserveDegrees()
	{
		List<String> genes = new ArrayList<String>(HPRD.getAllSymbols());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return ((Integer) HPRD.getInteractors(o2).size()).compareTo(HPRD.getInteractors(o1).size());
			}
		});

		int start = 0;
		int end = 0;
		int prevDegree = HPRD.getDegree(genes.get(0));

		List<String> shuffled = new ArrayList<String>();

		for (String gene : genes)
		{
			int degree = HPRD.getDegree(gene);

			if (degree == prevDegree)
			{
				end++;
			}
			else if (degree < prevDegree)
			{
				List<String> list = new ArrayList<String>(genes.subList(start, end));
				Collections.shuffle(list);
				shuffled.addAll(list);

				start = end;
				end++;
				prevDegree = degree;
			}
			else
			{
				throw new AssertionError("Should not reach here");
			}
		}
		List<String> list = genes.subList(start, end);
		Collections.shuffle(list);
		shuffled.addAll(list);

		assert genes.size() == shuffled.size();

		initShuffledMapping(genes, shuffled);
	}

	private void initShuffledAlmostPreserveDegrees()
	{
		List<String> genes = new ArrayList<String>(HPRD.getAllSymbols());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return ((Integer) HPRD.getInteractors(o2).size()).compareTo(HPRD.getInteractors(o1).size());
			}
		});

		int size = genes.size();
		int bins = 1000;
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

		initShuffledMapping(genes, shuffled);
	}

	private void initShuffledMapping(List<String> genes, List<String> shuffled)
	{
		redirectFwd = new HashMap<String, String>(genes.size());
		redirectBkw = new HashMap<String, String>(genes.size());

		for (int i = 0; i < genes.size(); i++)
		{
			redirectFwd.put(genes.get(i), shuffled.get(i));
			redirectBkw.put(shuffled.get(i), genes.get(i));
		}

		map = new HashMap<String, Set<String>>(genes.size());
	}

	private void writeDegrees(List<String> genes) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("/home/ozgun/Desktop/degrees.txt"));

		for (String gene : genes)
		{
			writer.write(HPRD.getInteractors(gene).size() + "\t" + gene + "\n");
		}

		writer.close();
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

		return new HashSet<String>(map.get(symbol));
	}
}
