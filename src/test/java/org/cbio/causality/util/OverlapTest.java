package org.cbio.causality.util;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Ozgun Babur
 */
public class OverlapTest
{
	@Test
	@Ignore
	public void testMutexAccuracy()
	{
		int n = 20;
		int a  = 10;
		int b = 11;

		System.out.println("overlap\tcalculated pval\tsimulated pval\tdifference");

		for (int o = Math.max(0, b-(n-a)); o <= Math.min(a, b); o++)
		{
			double e = (a * b) / (double) n;

			int trial = 100000;
			int hit = 0;

			for (int i = 0; i < trial; i++)
			{
				int ov = createRandomOverlap(n, a, b);
				if (ov <= o) hit++;
			}

			double sim = hit / (double) trial;
			double calc = Overlap.calcMutexPval(n, o, a, b);
			double dif = sim - calc;
			System.out.println(o + "\t" + sim + "\t" + calc + "\t" + dif);
			Assert.assertTrue(Math.abs(dif) < 0.01 || Math.abs(dif / calc) < 0.01);
		}
	}
	
	@Test
	@Ignore
	public void testMutexAccuracy2()
	{
		int n = 928;
		int a  = 305;
		int b = 284;
		int o = 67;

		System.out.println("overlap\tcalculated pval\tsimulated pval\tdifference");

		while((a+b-o) <= n)
		{
			int trial = 100000;
			int hit = 0;

			for (int i = 0; i < trial; i++)
			{
				int ov = createRandomOverlap(n, a, b);
				if (ov <= o) hit++;
			}

			double sim = hit / (double) trial;
			double calc = Overlap.calcMutexPval(n, o, a, b);
			double dif = sim - calc;
			System.out.println(o + "\t" + sim + "\t" + calc + "\t" + dif);
			Assert.assertTrue(Math.abs(dif) < 0.01 || Math.abs(dif / calc) < 0.01);
			a++;
			b++;
			o++;
		}
	}

	@Test
	@Ignore
	public void showTestValue()
	{
		int n = 158;
		int a  = 43;
		int b = 70;
		int o = 14;

		double pval = Overlap.calcMutexPval(n, o, a, b);
		System.out.println("pval = " + pval);
	}


	private static Random rand = new Random();

	private static int createRandomOverlap(int n, int a, int b)
	{
		List<Integer> list = new ArrayList<Integer>(n);
		for (int i = 0; i < n; i++)
		{
			list.add(i);
		}

		int overlap = 0;

		for (int i = 0; i < b; i++)
		{
			Integer x = list.remove(rand.nextInt(list.size()));
			if (x < a) overlap++;
		}

		return overlap;
	}
}
