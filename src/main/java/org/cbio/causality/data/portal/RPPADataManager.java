package org.cbio.causality.data.portal;

import org.cbio.causality.model.RPPAData;
import org.cbio.causality.util.ArrayUtil;

import java.io.IOException;
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

			RPPAData data = new RPPAData(s1[0], ArrayUtil.toDouble(s2),
				new HashSet<String>(Arrays.asList(s1[2].split("/"))),
				s1[1].equals("phosphorylation") ?
					parsePhosphoSites(s1[0], s1[3]) : Collections.<String>emptySet());

			for (String gene : data.genes)
			{
				if (!cache.containsKey(gene)) cache.put(gene, new ArrayList<RPPAData>());
				cache.get(gene).add(data);
			}
		}
	}

	private Set<String> parsePhosphoSites(String id, String pSite)
	{
		if (id.endsWith(pSite)) return Collections.singleton(pSite.substring(1));

		Set<String> sites = null;
		for (String s : id.split("_"))
		{
			if (s.contains("-")) s = s.substring(0, s.indexOf("-"));

			if (s.equals(pSite))
			{
				sites = new HashSet<String>();
				sites.add(pSite.substring(1));
			}
			else if (sites != null)
			{
				if (s.startsWith("p")) s = s.substring(1);
				sites.add(s);
			}
		}
		if (sites == null) return Collections.singleton(pSite.substring(1));
		return sites;
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
