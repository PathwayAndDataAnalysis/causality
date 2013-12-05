package org.cbio.causality.util;

import java.util.*;

/**
 * This class selects the subset of the given results to maintain a certain false discovery rate.
 * @author Ozgun Babur
 */
public class FDR
{
	public static List<String> select(final Map<String, Double> results, double fdrThr,
		List<Double> randomized, int randMultiplier)
	{
		Collections.sort(randomized);

		List<String> keys = new ArrayList<String>(results.keySet());

		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return results.get(o1).compareTo(results.get(o2));
			}
		});

		int ranIndex = 0;
		double ranPv = 0;
		int maxIndex = 0;

		for (int i = 0; i < keys.size(); i++)
		{
			String key = keys.get(i);
			double pval = results.get(key);

			while(ranPv <= pval && ranIndex < randomized.size())
			{
				ranPv = randomized.get(ranIndex++);
			}

			double noise = (ranIndex - 1) / randMultiplier;

			if (noise / (i+1) <= fdrThr) maxIndex = i;
		}

		return new ArrayList<String>(keys.subList(0, maxIndex+1));
	}

	/**
	 * @param results
	 * @param fdrThr
	 * @return
	 */
	public static List<String> select(final Map<String, Double> results, Map<String, Double> limits,
		double fdrThr)
	{
		List<Double> limitList = new ArrayList<Double>(limits.values());
		Collections.sort(limitList);

		List<String> keys = new ArrayList<String>(results.keySet());

		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return results.get(o1).compareTo(results.get(o2));
			}
		});

		int limIndex = 0;
		double limPv = limitList.get(0);
		int maxIndex = -1;

		for (int i = 0; i < keys.size(); i++)
		{
			String key = keys.get(i);
			double pval = results.get(key);

			while(limPv <= pval && limIndex < limitList.size()-1)
			{
				limPv = limitList.get(++limIndex);
			}

			double noise = pval * (limIndex + 1);

			if (noise / (i + 1) <= fdrThr) maxIndex = i;
		}

		if (maxIndex < 0) return Collections.emptyList();
		else return new ArrayList<String>(keys.subList(0, maxIndex+1));
	}

	/**
	 * @param results
	 * @return
	 */
	public static Map<String, Double> getQVals(final Map<String, Double> results,
		Map<String, Double> limits)
	{
		Map<String, Double> qvals = new HashMap<String, Double>();
		List<Double> limitList = new ArrayList<Double>(limits.values());
		Collections.sort(limitList);

		List<String> keys = new ArrayList<String>(results.keySet());

		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return results.get(o1).compareTo(results.get(o2));
			}
		});

		int limIndex = 0;
		double limPv = limitList.get(0);

		for (int i = 0; i < keys.size(); i++)
		{
			String key = keys.get(i);
			double pval = results.get(key);

			while(limPv <= pval && limIndex < limitList.size()-1)
			{
				limPv = limitList.get(++limIndex);
			}

			double noise = pval * (limIndex + 1);

			double fdr = noise / (i + 1);
			qvals.put(key, fdr);
		}
		return qvals;
	}
}
