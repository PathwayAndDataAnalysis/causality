package org.cbio.causality.util;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Consider below two-by-two table of counts of two traits.
 *
 *      A-  A+
 *     --------
 * B- | a | b |
 *    ---------
 * B+ | c | d |
 *    ---------
 *
 * Fisher's exact test calculates the exact probability of observing this pattern if traits A and B
 * were independent. This class provides two one-tailed tests, one for positive dependency and one
 * for negative dependency.
 *
 * @author Ozgun Babur
 */
public class FishersExactTest
{
	public static double calcNegativeDepPval(int a, int b, int c, int d)
	{
		double pval = 0;

		do {
			FactorialSolver f = new FactorialSolver(
				new ArrayList<Integer>(Arrays.asList(a+b, c+d, a+c, b+d)),
				new ArrayList<Integer>(Arrays.asList(a, b, c, d, a+b+c+d)));
			pval += f.solve();

			a--;
			b++;
			c++;
			d--;
		}
		while(d >= 0 && a >= 0);
		return pval;
	}

	public static double calcPositiveDepPval(int a, int b, int c, int d)
	{
		double pval = 0;

		do {
			FactorialSolver f = new FactorialSolver(
				new ArrayList<Integer>(Arrays.asList(a+b, c+d, a+c, b+d)),
				new ArrayList<Integer>(Arrays.asList(a, b, c, d, a+b+c+d)));
			pval += f.solve();

			a++;
			b--;
			c--;
			d++;
		}
		while(b >= 0 && c >= 0);
		return pval;
	}

	public static double calcEnrichmentPval(int size, int featuredOverall, int selected,
		int featuredSelected)
	{
		assert selected <= size;
		assert featuredSelected <= selected;
		assert featuredSelected <= featuredOverall;

		return calcPositiveDepPval(size - selected - featuredOverall + featuredSelected,
			featuredOverall - featuredSelected, selected - featuredSelected, featuredSelected);
	}

	public static void main(String[] args)
	{
	}
}
