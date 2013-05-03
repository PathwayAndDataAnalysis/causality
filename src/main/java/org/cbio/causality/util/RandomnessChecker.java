package org.cbio.causality.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * This class is useful when we want to check is a set of p-values belong to random events.
 * @author Ozgun Babur
 */
public class RandomnessChecker
{
	/**
	 * Set of p-values.
	 */
	List<Double> values;

	/**
	 * Constructor.
	 */
	public RandomnessChecker()
	{
		this(1000);
	}

	/**
	 * Constructor.
	 * @param initialSize initial size of the values list
	 */
	public RandomnessChecker(int initialSize)
	{
		values = new ArrayList<Double>(initialSize);
	}

	/**
	 * Adds another p-value to the p-value set.
	 * @param pval p-value to add
	 */
	public void add(double pval)
	{
		values.add(pval);
	}

	/**
	 * Gets the number of p-values in the set that are smaller than the given threshold.
	 * @param thr threshold value
	 * @return number of smaller p-values
	 */
	public int getCount(double thr)
	{
		int cnt = 0;
		for (Double value : values)
		{
			if (value <= thr) cnt++;
		}
		return cnt;
	}

	/**
	 * Prepares a status string for the given threshold. Shows expected, observed, difference, and
	 * diff percentage values.
	 * @param thr threshold value
	 * @return summary string
	 */
	public String getStatusForThreshold(double thr)
	{
		double e = getExpected(thr);
		int o = getCount(thr);
		double d = o - e;
		return "Expected: " + e + "   Observed: " + o + "   Diff: " + (d > 0 ? "+" : "") + d +
			"   Off by " + (int)(Math.abs(d) * 100 / e) + "%";
	}

	/**
	 * Gets the expected number of p-values in the set that are smaller than the given threshold,
	 * assuming that the p-values are randomly drawn.
	 * @param thr threshold value
	 * @return expected number of p-values
	 */
	public double getExpected(double thr)
	{
		return values.size() * thr;
	}

	public static void main(String[] args)
	{
		RandomnessChecker rc = new RandomnessChecker();
		for (int i = 0; i < 10000; i++)
		{
			rc.add(Math.random());
		}
		System.out.println(rc.getStatusForThreshold(0.05));
	}
}
