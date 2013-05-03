package org.cbio.causality.util;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

/**
 * Performs chi-square related operations.
 * @author Ozgun Babur
 */
public class ChiSquare
{
	/**
	 * Calculates the p-value of the given chi-square value with the given degrees of freedom.
	 * @param x chi value
	 * @param n degrees of freedom
	 * @return p-value
	 */
	public static double pValue(double x, double n)
	{
		if(n==1 && x>1000)
		{
			return 0;
		}
		if(x>1000 || n>1000)
		{
			double q = pValue((x - n) * (x - n) / (2 * n), 1) / 2;

			if(x>n)
			{
				return q;
			}
			else
			{
				return 1-q;
			}
		}
		double p = Math.exp(-0.5 * x);
		if((n % 2) == 1)
		{
			p = p * Math.sqrt(2 * x / Math.PI);
		}

		double k = n;

		while(k >= 2)
		{
			p = p * x / k;
			k = k - 2;
		}
		double t = p;
		double a = n;
		while(t > 0.0000000001 * p)
		{
			a = a + 2;
			t = t * x / a;
			p = p + t;
		}

		return 1 - p;
	}

	public static void main(String[] args) throws MathException
	{
		RandomnessChecker rc = new RandomnessChecker();
		for (int i = 0; i < 10000; i++)
		{
			rc.add(Math.random());
		}
		System.out.println(rc.getStatusForThreshold(0.05));
	}
}
