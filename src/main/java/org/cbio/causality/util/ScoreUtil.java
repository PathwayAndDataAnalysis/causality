package org.cbio.causality.util;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class ScoreUtil
{
	private List<Double> scores;
	private double max;
	private double min;
	protected int multiFactor;

	public ScoreUtil()
	{
		this(1);
	}

	public ScoreUtil(int multiFactor)
	{
		this.multiFactor = multiFactor;
		scores = new ArrayList<Double>();
	}

	public void addSCore(double d)
	{
		if (scores.isEmpty())
		{
			min = d;
			max = d;
		}

		scores.add(d);

		if (d > max) max = d;
		else if (d < min) min = d;
	}

	public void unite(ScoreUtil sc)
	{
		this.scores.addAll(sc.scores);
		this.multiFactor += sc.multiFactor;
	}

	public int countOverThr(double thr)
	{
		int cnt = 0;
		for (Double score : scores)
		{
			if (score >= thr) cnt++;
		}
		return cnt;
	}

	public static final double EPS = 0.00000001;

	/**
	 * Performs a logarithmic search for the given FDR.
	 * @param real scores for the non-random case(s)
	 * @param fdr target FDR
	 * @return threshold value for the target FDR
	 */
	public double getThresholdForFDR(ScoreUtil real, double fdr)
	{
		if (getFDRForThr(real, real.max) > fdr) return real.max + EPS;

		double highFDR = getFDRForThr(real, this.min);
		if (highFDR <= fdr) return this.min;

		double max = this.max;
		double min = this.min;

		if (min == max) return min + EPS;

		double lowFDR = getFDRForThr(real, this.max + EPS);
		double prevValForHighFDR = this.min;
		double prevValForLowFDR = this.max + EPS;

		double val = (min + max) / 2;
		double FDR = getFDRForThr(real, val);

		do
		{
			double temp = val;
			if (FDR >= fdr)
			{
				val = (val + max) / 2;
				min = temp;
			}
			else
			{
				val = (val + min) / 2;
				max = temp;
			}

			if (FDR > fdr)
			{
				highFDR = FDR;
				prevValForHighFDR = temp;
			}
			else if (FDR < fdr)
			{
				lowFDR = FDR;
				prevValForLowFDR = temp;
			}

			FDR = getFDRForThr(real, val);
			System.out.println("FDR = " + FDR);
		}
		while (highFDR >= fdr && lowFDR < fdr && Math.abs(prevValForLowFDR - prevValForHighFDR) < EPS);

		return val;
	}

	public double getFDRForThr(ScoreUtil real, double thr)
	{
		int ranCnt = countOverThr(thr);
		int norCnt = real.countOverThr(thr);

		if (norCnt == 0) return 0;

		double ranNormalized = ranCnt / (double) multiFactor;
		double norNormalized = norCnt / (double) real.multiFactor;

		return ranNormalized / norNormalized;
	}
}
