package org.cbio.causality.data.portal;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Reads and caches expression data from cBioPortal.
 *
 * @author Ozgun Babur
 */
public class ExpDataManager
{
	private GeneticProfile profile;
	private CaseList caseList;

	private Map<String, double[]> cache;
	private CBioPortalManager cman;
	private Set<String> notFound;

	public ExpDataManager(GeneticProfile profile, CaseList caseList)
	{
		this.profile = profile;
		this.caseList = caseList;

		cache = new HashMap<String, double[]>();
		cman = new CBioPortalManager();
		notFound = new HashSet<String>();
	}

	public double[] get(String symbol)
	{
		if (notFound.contains(symbol)) return null;

		if (!cache.containsKey(symbol))
		{
			String[] data = cman.getDataForGene(symbol, profile, caseList);

			if (data != null) cache.put(symbol, stringToDouble(data));
		}

		if (cache.containsKey(symbol))
		{
			return cache.get(symbol);
		}
		else
		{
			notFound.add(symbol);
			return null;
		}
	}

	private double[] stringToDouble(String[] data)
	{
		double[] d = new double[data.length];

		for (int i = 0; i < data.length; i++)
		{
			try
			{
				d[i] = log2(Double.parseDouble(data[i]));
//				d[i] = Double.parseDouble(data[i]);
			}
			catch (NumberFormatException e)
			{
				d[i] = Double.NaN;
			}
		}
		return d;
	}

	public static final double LOG2 = Math.log(2);
	private double log2(double val)
	{
		return Math.log(val) / LOG2;
	}
}
