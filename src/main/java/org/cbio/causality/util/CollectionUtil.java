package org.cbio.causality.util;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class CollectionUtil
{
	public static List<String> getSortedList(final Map<String, Double> map)
	{
		List<String> list = new ArrayList<String>(map.keySet());
		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				Double v1 = map.get(o1);
				Double v2 = map.get(o2);

				if (v1.equals(v2)) return o1.compareTo(o2);
				return v1.compareTo(v2);
			}
		});

		return list;
	}
}
