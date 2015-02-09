package org.cbio.causality.util;

import org.apache.commons.math.MathException;

/**
 * Aggregates independent p-values using Fisher's combined probability test.
 * @author Ozgun Babur
 */
public class FishersCombinedProbability
{
	/**
	 * Calculates the combined p-value for the given independent p-values.
	 * @param pvals p-values to combine
	 * @return aggregate p-value
	 */
	public static double pValue(double... pvals)
	{
		double chi = 0;

		for (double pval : pvals)
		{
			chi += -2 * Math.log(pval);
		}

		return ChiSquare.pValue(chi, 2 * pvals.length);
	}

	public static void main(String[] args) throws MathException
	{
		System.out.println(pValue(0.9, 0.9, 0.9, 0.9, 0.1));
	}
}
