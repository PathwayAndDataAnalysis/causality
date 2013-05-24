package org.cbio.causality.hprd;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ShuffledHPRD implements InteractionProvider
{
	private Map<String, String> redirectFwd;
	private Map<String, String> redirectBkw;
	private Map<String, Set<String>> map;

	public ShuffledHPRD()
	{
		List<String> list1 = new ArrayList<String>(HPRD.getAllSymbols());
		List<String> list2 = new ArrayList<String>(list1);
		Collections.shuffle(list2);

		redirectFwd = new HashMap<String, String>(list1.size());
		redirectBkw = new HashMap<String, String>(list1.size());

		for (int i = 0; i < list1.size(); i++)
		{
			redirectFwd.put(list1.get(i), list2.get(i));
			redirectBkw.put(list2.get(i), list1.get(i));
		}

		map = new HashMap<String, Set<String>>(list1.size());
	}

	public Set<String> getInteractions(String symbol)
	{
		if (map.containsKey(symbol)) return map.get(symbol);

		String s = redirectBkw.get(symbol);

		Set<String> inter = HPRD.getInteractors(s);

		Set<String> redirInter = new HashSet<String>(inter.size());

		for (String s1 : inter)
		{
			redirInter.add(redirectFwd.get(s1));
		}
		map.put(symbol, redirInter);

		return redirInter;
	}
}
