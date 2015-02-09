package org.cbio.causality.util;

import org.cbio.causality.util.trendline.PolyTrendLine;

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
		if (results.isEmpty()) return Collections.emptyList();
		if (randomized.isEmpty()) return new ArrayList<String>(results.keySet());

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
		int maxIndex = -1;

		for (int i = 0; i < keys.size(); i++)
		{
			String key = keys.get(i);
			double pval = results.get(key);

			while(ranPv <= pval && ranIndex < randomized.size())
			{
				ranPv = randomized.get(ranIndex++);
			}

			double noise = (ranIndex - 1) / (double) randMultiplier;

			if (noise / (i+1) <= fdrThr) maxIndex = i;
		}

		if (maxIndex < 0) return Collections.emptyList();
		else return new ArrayList<String>(keys.subList(0, maxIndex+1));
	}

	/**
	 * @param results
	 * @param fdrThr
	 * @return
	 */
	public static List<String> select(Map<String, Double> results, Map<String, Double> limits,
		final Map<String, Double> priorityScores, double fdrThr)
	{
		List<String> items = new ArrayList<String>(priorityScores.keySet());
		Collections.sort(items, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return priorityScores.get(o2).compareTo(priorityScores.get(o1));
			}
		});

		if (limits == null) limits = defaultLimits(results);

		int priorityIndex = -1;
		int maxResult = 0;
		int selectedPriority = -1;
		List<String> result = null;

		do
		{
			while(priorityIndex < items.size() - 1 && (priorityIndex < 0 ||
				priorityScores.get(items.get(priorityIndex)).equals(
					priorityScores.get(items.get(priorityIndex + 1)))))
			{
				priorityIndex++;
			}

			Map<String, Double> resultsTemp = new HashMap<String, Double>();
			Map<String, Double> limitsTemp = new HashMap<String, Double>();

			for (int i = 0; i <= priorityIndex; i++)
			{
				String key = items.get(i);
				resultsTemp.put(key, results.get(key));
				limitsTemp.put(key, limits.get(key));
			}

			List<String> selected = select(resultsTemp, limitsTemp, fdrThr);

			if (selected.size() >= maxResult)
			{
				maxResult = selected.size();
				selectedPriority = priorityIndex;
				result = selected;
			}

			priorityIndex++;
		}
		while(priorityIndex < items.size() - 1);

		System.out.println("selectedPriority = " + selectedPriority);

		if (result == null) return Collections.emptyList();
		return result;
	}

	public static List<String> selectBH(final Map<String, Double> results, double fdrThr)
	{
		return select(results, null, fdrThr);
	}

	public static List<String> select(final Map<String, Double> results, Map<String, Double> limits,
		double fdrThr)
	{
		if (results.isEmpty()) return Collections.emptyList();

		if (limits == null) limits = defaultLimits(results);

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

	private static Map<String, Double> defaultLimits(Map<String, Double> results)
	{
		Map<String, Double> limits;
		limits = new HashMap<String, Double>(results);

		for (String key : new HashSet<String>(limits.keySet()))
		{
			limits.put(key, 0D);
		}
		return limits;
	}

	public static List<String> selectWithPvalThreshold(final Map<String, Double> pvals, double pvalThr)
	{
		List<String> keys = new ArrayList<String>(pvals.keySet());
		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Double(pvals.get(o1)).compareTo(new Double(pvals.get(o2)));
			}
		});

		int cut = 0;
		while (pvals.get(keys.get(cut)) <= pvalThr && cut < keys.size()) cut++;

		return keys.subList(0, cut);
	}

	/**
	 * @param results
	 * @return
	 */
	public static Map<String, Double> getQVals(final Map<String, Double> results,
		Map<String, Double> limits)
	{
		Map<String, Double> qvals = new HashMap<String, Double>();

		if (limits == null) limits = defaultLimits(results);

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

	public static int[] getResultSizesUsingPolyCurve(final Map<String, Double> results, double[] thrs)
	{
		int[] sizes = new int[thrs.length];
		for (int i = 0; i < thrs.length; i++)
		{
			sizes[i] = selectUsingPolyCurve(results, thrs[i]).size();
		}
		return sizes;
	}

	/**
	 * @param results
	 * @param fdrThr
	 * @return
	 */
	public static List<String> selectUsingPolyCurve(final Map<String, Double> results, double fdrThr)
	{
		if (results.isEmpty()) return Collections.emptyList();

		List<String> keys = new ArrayList<String>(results.keySet());

		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return results.get(o1).compareTo(results.get(o2));
			}
		});

		double noiseSize = estimateNoiseVolume(results);

		int maxIndex = -1;

		for (int i = 0; i < keys.size(); i++)
		{
			String key = keys.get(i);
			double pval = results.get(key);

			double noise = pval * noiseSize;

			if (noise / (i + 1) <= fdrThr) maxIndex = i;
		}

		if (maxIndex < 0) return Collections.emptyList();
		else return new ArrayList<String>(keys.subList(0, maxIndex+1));
	}

	public static double estimateNoiseVolume(Map<String, Double> pvals)
	{
		double[][] f = getNoiseEstimatesForDifferentLambda(pvals);

		PolyTrendLine trendLine = new PolyTrendLine(3);
		trendLine.setValues(f[1], f[0]);

		double noiseRatio = trendLine.predict(f[0][f[0].length - 1]);

//		System.out.println("\n\nplot");
//		for (int i = 0; i < f[0].length; i++)
//		{
//			System.out.println(f[0][i] + "\t" + f[1][i] + "\t" + trendLine.predict(f[0][i]));
//		}
//		System.out.println();

		return noiseRatio * pvals.size();
	}

	public static double[][] getNoiseEstimatesForDifferentLambda(Map<String, Double> pvals)
	{
		List<Double> vals = new ArrayList<Double>();

		int[] cnt = new int[99];

		for (Double val : pvals.values())
		{
			for (int i = 0; i < cnt.length; i++)
			{
				if (val > ((i + 1) * 0.01)) cnt[i]++;
			}
		}

		double m = pvals.size();

		for (int i = 0; i < cnt.length; i++)
		{
			double x = cnt[i] / (m * (1D - ((i + 1) * 0.01)));

			if (x > 0) vals.add(x);
			else break;
		}

		double[][] v = new double[2][vals.size()];

		for (int i = 0; i < v[0].length; i++)
		{
			v[0][i] = (i + 1) * 0.01;
			v[1][i] = vals.get(i);
		}

		return v;
	}

	public static double decideBestFDR(final  Map<String,  Double>  results,
			List<Double>  randomized, int randMultiplier)
	{
		double bestFDR = -1;
		double maxScore = 0;

		System.out.println("\nFDR\tResult size\tExpected true positives\ttp-fp");
		for (int i = 1; i <= 50; i++)
		{
			double fdr = i / 100D;
			List<String> select = select(results, fdr, randomized, randMultiplier);
			double tp = select.size() * (1 - fdr);
			double fp = select.size() * fdr;
			double score = tp - fp;
			if (score > maxScore)
			{
				maxScore = score;
				bestFDR = fdr;
			}

			System.out.println(fdr + "\t" + select.size() + "\t" + ((int) Math.round(tp)) + "\t" + ((int) Math.round(tp - fp)));
		}
		System.out.println();

		return bestFDR;
	}

	public static double decideBestFDR_BH(final  Map<String,  Double>  results, Map<String,  Double>  limits)
	{
		double bestFDR = -1;
		double maxScore = 0;

//		System.out.println("\nFDR\tResult size\tExpected true positives\ttp-fp");
		for (int i = 1; i <= 50; i++)
		{
			double fdr = i / 100D;
			List<String> select = select(results, limits, fdr);
			double tp = select.size() * (1 - fdr);
			double fp = select.size() * fdr;
			double score = tp - fp;
			if (score > maxScore)
			{
				maxScore = score;
				bestFDR = fdr;
			}

//			System.out.println(fdr + "\t" + select.size() + "\t" + ((int) Math.round(tp)) + "\t" + ((int) Math.round(tp - fp)));
		}
//		System.out.println();

		return bestFDR;
	}
}
