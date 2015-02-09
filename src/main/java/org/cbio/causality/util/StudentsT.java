package org.cbio.causality.util;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.TestUtils;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class StudentsT
{
	public static double getPValOfMeanDifference(double[] x0, double[] x1)
	{
		if (x0.length < 2 || x1.length < 2) return Double.NaN;

		try
		{
			return TestUtils.tTest(x0, x1);
		}
		catch (MathException e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	public static double getPValOfMeanEqualTo(double[] x, double mean)
	{
		if (x.length < 2) return Double.NaN;

		try
		{
			return TestUtils.tTest(mean, x);
		}
		catch (MathException e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	public static double getPValOfMeanDifference_OldAndWrong(double[] x0, double[] x1)
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
		Histogram h = new Histogram(0.01);
		h.setBorderAtZero(true);
		for (int j = 0; j <200000; j++)
		{
			List<Double> list1 = new ArrayList<Double>();
			List<Double> list2 = new ArrayList<Double>();

			Random rand = new Random();

			for (int i = 0; i < 20; i++)
			{
				list1.add(rand.nextGaussian());
			}
			for (int i = 0; i < 5; i++)
			{
				list2.add(rand.nextGaussian());
			}

			double[] vals1 = ArrayUtil.toArray(list1, 0D);
			double[] vals2 = ArrayUtil.toArray(list2, 0D);

			double p = getPValOfMeanDifference(vals1, vals2);

			h.count(p);
		}

		h.print();
	}

	public static double getPValOfMeanDifferenceBySimulation(double[] x0, double[] x1,
		int maxTrials, int maxHits)
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

		int trial;
		for (trial = 0; trial < maxTrials; trial++)
		{
			Collections.shuffle(nums);

			double d = getDifOfMean(nums, xx0, xx1);
			if (d >= dif) hit++;

			if (hit == maxHits)
			{
				trial++;
				break;
			}
		}

		return hit / (double) trial;
	}

//	public static double[] getPValOfMeanDifferenceBySimulation(double[] x0, double[] x1,
//		int maxTrials, int maxHits)
//	{
//		double mean0 = Summary.mean(x0);
//		double mean1 = Summary.mean(x1);
//
//		double dif = Math.abs(mean0 - mean1);
//
//		List<Double> nums = new ArrayList<Double>(x0.length + x1.length);
//		for (double v : x0) nums.add(v);
//		for (double v : x1) nums.add(v);
//
//		int hit = 0;
//		double[] xx0 = new double[x0.length];
//		double[] xx1 = new double[x1.length];
//
//		Collections.sort(nums);
//		double maxDif = getDifOfMean(nums, xx0, xx1);
//		int hitOfMax = 0;
//
//		int trial;
//		for (trial = 0; trial < maxTrials; trial++)
//		{
//			Collections.shuffle(nums);
//
//			double d = getDifOfMean(nums, xx0, xx1);
//			if (d >= dif) hit++;
//
//			if (d == maxDif) hitOfMax++;
//			if (hit == maxHits)
//			{
//				trial++;
//				break;
//			}
//		}
//
//		return new double[]{hit / (double) trial, hitOfMax / (double) trial};
//	}

	private static double getDifOfMean(List<Double> nums, double[] xx0, double[] xx1)
	{
		for (int j = 0; j < xx0.length; j++)
		{
			xx0[j] = nums.get(j);
		}
		for (int j = 0; j < xx1.length; j++)
		{
			xx1[j] = nums.get(xx0.length + j);
		}

		return Math.abs(Summary.mean(xx0) - Summary.mean(xx1));
	}
}
