package org.cbio.causality.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Common functions related to arrays.
 *
 * @author Ozgun Babur
 */
public class ArrayUtil
{
	public static boolean[] negate(boolean[] posit)
	{
		boolean[] neg = new boolean[posit.length];
		for (int i = 0; i < posit.length; i++)
		{
			neg[i] = !posit[i];
		}
		return neg;
	}

	public static int countValue(boolean[] b, boolean val)
	{
		int cnt = 0;
		for (boolean v : b)
		{
			if (v == val) cnt++;
		}
		return cnt;
	}

	/**
	 * Count of indexes where both arrays have the value.
	 */
	public static int countValue(boolean[] b1, boolean[] b2, boolean val)
	{
		int cnt = 0;
		for (int i = 0; i < b1.length; i++)
		{
			if (b1[i] == val && b2[i] == val) cnt++;
		}
		return cnt;
	}

	public static void ORWith(boolean[] toChange, boolean[] toAdd)
	{
		if (toChange.length != toAdd.length) throw new IllegalArgumentException(
			"Array sizes have to be equal.");

		for (int i = 0; i < toAdd.length; i++)
		{
			if (toAdd[i]) toChange[i] = true;
		}
	}

	public static double[] subset(double[] vals, boolean[] select)
	{
		List<Double> list = new ArrayList<Double>();
		for (int i = 0; i < vals.length; i++)
		{
			if (select[i] && !Double.isNaN(vals[i])) list.add(vals[i]);
		}
		double[] sub = new double[list.size()];

		int i = 0;
		for (Double val : list)
		{
			sub[i++] = val;
		}
		return sub;
	}

	public static <T> boolean[] getLocations(T[] array, T query)
	{
		boolean[] loc = new boolean[array.length];
		for (int i = 0; i < array.length; i++)
		{
			loc[i] = array[i].equals(query);
		}
		return loc;
	}

	public static int[] toArray(List<Integer> vals, int dummy)
	{
		int[] array = new int[vals.size()];
		int i = 0;
		for (Integer val : vals)
		{
			array[i++] = val;
		}
		return array;
	}

	public static double[] toArray(List<Double> vals, double dummy)
	{
		double[] array = new double[vals.size()];
		int i = 0;
		for (Double val : vals)
		{
			array[i++] = val;
		}
		return array;
	}

	public static double[] toDouble(String[] s)
	{
		double[] v = new double[s.length];
		for (int i = 0; i < s.length; i++)
		{
			try
			{
				v[i] = Double.parseDouble(s[i]);
			}
			catch (NumberFormatException e)
			{
				v[i] = Double.NaN;
			}
		}
		return v;
	}

	public static double[] toPrimitive(Double[] vals)
	{
		double[] v = new double[vals.length];
		for (int i = 0; i < v.length; i++)
		{
			v[i] = vals[i];
		}
		return v;
	}
}
