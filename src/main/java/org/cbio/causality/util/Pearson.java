package org.cbio.causality.util;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;

import java.util.Random;

/**
 * @author Ozgun Babur
 */
public class Pearson
{
	public static double correlation(double[] v1, double[] v2)
	{
		PearsonsCorrelation pc = new PearsonsCorrelation();
		return pc.correlation(v1, v2);
	}

	public static double corrPval(double[] v1, double[] v2) throws MathException
	{
		PearsonsCorrelation pc = new PearsonsCorrelation(transform(v1, v2));
		return pc.getCorrelationPValues().getColumn(0)[1];
	}

	public static double[][] transform(double[] v1, double[] v2)
	{
		double[][] d = new double[v1.length][2];
		for (int i = 0; i < v1.length; i++)
		{
			d[i][0] = v1[i];
			d[i][1] = v2[i];
		}
		return d;
	}

	public static double spearmanCorrelation(double[] v1, double[] v2)
	{
		SpearmansCorrelation pc = new SpearmansCorrelation();
		return pc.correlation(v1, v2);
	}

	public static void main(String[] args) throws MathException
	{
		Random r = new Random();
		int size = 20;
		double[] v1 = new double[size];
		double[] v2 = new double[size];
		for (int i = 0; i < size; i++)
		{
			v1[i] = r.nextDouble();
			v2[i] = r.nextDouble();
		}

		System.out.println("correlation(v1, v2) = " + correlation(v1, v2));
		System.out.println("corrPval(v1, v2) = " + corrPval(v1, v2));
	}
}
