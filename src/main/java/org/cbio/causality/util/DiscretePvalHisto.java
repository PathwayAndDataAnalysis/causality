package org.cbio.causality.util;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class DiscretePvalHisto
{
	double[] pvals;
	double window;

	public DiscretePvalHisto(double[] pvals, double window)
	{
		this.pvals = pvals;
		this.window = window;
		retouch();
	}

	public DiscretePvalHisto(List<Double> pvalList, double window)
	{
		this.pvals = new double[pvalList.size()];
		for (int i = 0; i < pvals.length; i++)
		{
			pvals[i] = pvalList.get(i);
		}
		this.window = window;
		retouch();
	}

	private void retouch()
	{
		for (int i = 0; i < pvals.length; i++)
		{
			pvals[i] = Math.round(pvals[i] * 1E10) / 1E10;
		}
	}

	public void plot()
	{
		Map<Double, Integer> map = new HashMap<Double, Integer>();

		for (double v : pvals)
		{
			if (!map.containsKey(v)) map.put(v, 0);
			map.put(v, map.get(v) + 1);
		}

		int total = 0;
		for (Integer cnt : map.values())
		{
			total += cnt;
		}

		List<Double> list = new ArrayList<Double>(map.keySet());
		Collections.sort(list);

		int[] border = getBorderIndexes(list);

		int from = -1;
		for (int b : border)
		{
			double cnt = 0;
			for (int j = from + 1; j <= b; j++)
			{
				cnt += map.get(list.get(j));
			}

			double dist = list.get(b) - (from >= 0 ? list.get(from) : 0);

			double density = cnt / (dist * total);

			System.out.println(list.get(b) + "\t" + density);
			from = b;
		}
	}

	private int[] getBorderIndexes(List<Double> list)
	{
		LinkedList<Integer> bord = new LinkedList<Integer>();

		double target = window;
		for (int i = 0; i < list.size(); i++)
		{
			if (list.get(i) < target) continue;

			bord.add(i);

			while (target < list.get(bord.getLast()) + window) target += window;
		}

		if (!bord.isEmpty() && bord.getLast() < list.size() - 1) bord.add(list.size() - 1);

		int i = 0;
		int[] inds = new int[bord.size()];
		for (Integer integ : bord.toArray(new Integer[bord.size()]))
		{
			inds[i++] = integ;
		}
		return inds;
	}
}
