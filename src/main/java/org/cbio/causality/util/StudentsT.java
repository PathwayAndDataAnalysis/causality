package org.cbio.causality.util;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class StudentsT
{
	public static double getPValOfMeanDifference(double[] x0, double[] x1)
	{
		if (x0.length == 0 || x1.length == 0) return 1;

		double mean0 = Summary.mean(x0);
		double mean1 = Summary.mean(x1);

		double var0 = Summary.variance(x0);
		double var1 = Summary.variance(x1);

		if (Double.isNaN(mean0) || Double.isNaN(mean1)) return 1;

		double sd = calcSDForEqualVar(x0.length, x1.length, var0, var1);

		if (sd == 0) return 1;

		double v = Math.abs(mean0 - mean1) / (sd * SQRT2);
		return 1 - ErrorFunction.getSignif(v);
	}

	private static double calcSDForEqualVar(int n0, int n1, double var0, double var1)
	{
		double sd = (var0 * (n0 - 1)) + (var1 * (n1 - 1));
		sd /= n0 + n1 - 2;
		sd = Math.sqrt(sd);
		sd *= Math.sqrt((1D / n0) + (1D / n1));
		return sd;
	}

	private static double calcSDForUnequalVar(int n0, int n1, double var0, double var1)
	{
		double sd = (var0 / n0) + (var1 / n1);
		sd = Math.sqrt(sd);
		return sd;
	}

	/**
	 * Square root of 2, to use in calculations.
	 */
	public static final double SQRT2 = Math.sqrt(2);

	public static void main(String[] args)
	{
		int iter = 10000;
		double[] p = new double[iter];
		double[] s = new double[iter];
		for (int k = 0; k < p.length; k++)
		{
			Random r = new Random();
			double[] x0 = new double[3];
			double[] x1 = new double[3];

			for (int i = 0; i < x0.length; i++)
			{
				x0[i] = r.nextGaussian();
			}
			for (int i = 0; i < x1.length; i++)
			{
				x1[i] = r.nextGaussian();
			}
			p[k] = getPValOfMeanDifference(x0, x1);
			s[k] = getPValOfMeanDifferenceBySimulation(x0, x1, 1000);
		}

		int cnt1 = 0;
		for (double v : p)
		{
			if (v < 0.05) cnt1++;
		}
		int cnt2 = 0;
		for (double v : s)
		{
			if (v < 0.05) cnt2++;
		}
		System.out.println("ratio = " + cnt1/(double)p.length);
		System.out.println("ratio = " + cnt2/(double)s.length);
	}

	public static double getPValOfMeanDifferenceBySimulation(double[] x0, double[] x1, int trials)
	{
		double mean0 = Summary.mean(x0);
		double mean1 = Summary.mean(x1);

		double dif = Math.abs(mean0 - mean1);

		List<Double> nums = new ArrayList<Double>(x0.length + x1.length);
		for (double v : x0) nums.add(v);
		for (double v : x1) nums.add(v);

		int hit = 0;
		double[] xx0 = new double[x0.length];
		double[] xx1 = new double[x1.length];

		for (int i = 0; i < trials; i++)
		{
			Collections.shuffle(nums);

			for (int j = 0; j < x0.length; j++)
			{
				xx0[j] = nums.get(j);
			}
			for (int j = 0; j < x1.length; j++)
			{
				xx1[j] = nums.get(x0.length + j);
			}

			double d = Math.abs(Summary.mean(xx0) - Summary.mean(xx1));
			if (d >= dif) hit++;
		}

		return hit / (double) trials;
	}
}
