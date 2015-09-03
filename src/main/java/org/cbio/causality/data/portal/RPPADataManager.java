package org.cbio.causality.data.portal;

import org.cbio.causality.rppa.RPPAData;
import org.cbio.causality.util.ArrayUtil;

import java.util.*;

/**
 * Reads and caches RPPA data from cBioPortal.
 *
 * @author Ozgun Babur
 */
public class RPPADataManager
{
	private CaseList caseList;

	private Map<String, List<RPPAData>> cache;
	private CBioPortalManager cman;

	public RPPADataManager(CaseList caseList)
	{
		this.caseList = caseList;

		cache = new HashMap<String, List<RPPAData>>();
		cman = new CBioPortalManager();
		load();
	}

	private void load()
	{
		Map<String[], String[]> map = cman.getRPPAData(caseList);

		for (String[] s1 : map.keySet())
		{
			String[] s2 = map.get(s1);

			List<String> genes = Arrays.asList(s1[2].split("/"));
			RPPAData data = new RPPAData(s1[0], new double[][]{ArrayUtil.toDouble(s2)},
				genes, s1[1].equals("phosphorylation") ?
					parsePhosphoSites(s1[0], genes, s1[3]) : null);

			for (String gene : data.genes)
			{
				if (!cache.containsKey(gene)) cache.put(gene, new ArrayList<RPPAData>());
				cache.get(gene).add(data);
			}
		}
	}

	private Map<String, List<String>> parsePhosphoSites(String id, List<String> genes, String pSite)
	{
		List<String> list = new ArrayList<String>();
		if (id.endsWith(pSite)) list.add(pSite.substring(1));
		else
		{
			List<String> sites = null;
			for (String s : id.split("_"))
			{
				if (s.contains("-")) s = s.substring(0, s.indexOf("-"));

				if (s.equals(pSite))
				{
					sites = new ArrayList<String>();
					sites.add(pSite.substring(1));
				}
				else if (sites != null)
				{
					if (s.startsWith("p")) s = s.substring(1);
					sites.add(s);
				}
			}
			if (sites == null) list.add(pSite.substring(1));
			else list.addAll(sites);
		}

		Map<String, List<String>> map = new HashMap<String, List<String>>();
		for (String gene : genes)
		{
			map.put(gene, list);
		}
		return map;
	}

	public List<RPPAData> get(String gene)
	{
		return cache.get(gene);
	}

	public boolean contains(String gene)
	{
		return cache.containsKey(gene);
	}

	public Set<String> getSymbols()
	{
		return cache.keySet();
	}

	public Set<RPPAData> getAllData()
	{
		Set<RPPAData> data = new HashSet<RPPAData>();
		for (List<RPPAData> list : cache.values())
		{
			data.addAll(list);
		}
		return data;
	}
}
