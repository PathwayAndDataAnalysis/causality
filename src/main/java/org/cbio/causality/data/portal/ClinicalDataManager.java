package org.cbio.causality.data.portal;

import org.cbio.causality.util.ArrayUtil;

import java.util.*;

/**
 * Reads and caches expression data from cBioPortal.
 *
 * @author Ozgun Babur
 */
public class ClinicalDataManager
{
	private Map<String, Double> os;
	private Map<String, Double> dfs;

	public ClinicalDataManager(CBioPortalAccessor acc)
	{
		this(acc.getCurrentCaseList());
	}

	public ClinicalDataManager(CaseList caseList)
	{
		CBioPortalManager cman = new CBioPortalManager();
		Map<String, Double>[] data = cman.getClinicalData(caseList);
		if (data != null)
		{
			dfs = data[0];
			os = data[1];
		}
	}

	public double[] getDFS(Set<String> caseIDs)
	{
		return get(caseIDs, dfs);
	}

	public double[] getOS(Set<String> caseIDs)
	{
		return get(caseIDs, os);
	}

	private double[] get(Set<String> caseIDs, Map<String, Double> map)
	{
		List<Double> vals = new ArrayList<Double>();
		for (String caseID : caseIDs)
		{
			Double val = map.get(caseID);
			if (val != null && !Double.isNaN(val)) vals.add(val);
		}
		return ArrayUtil.toArray(vals, 0);
	}
}
