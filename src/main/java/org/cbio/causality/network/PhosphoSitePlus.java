package org.cbio.causality.network;

import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.TermCounter;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PhosphoSitePlus
{
	static Map<String, Map<String, Integer>> typeMap;
	static Map<String, Map<String, String>> actualMap;

	public static Integer getEffect(String gene, String site)
	{
		if (typeMap.containsKey(gene))
		{
			return typeMap.get(gene).get(site);
		}
		return null;
	}

	public static Integer getClosestEffect(String gene, String site)
	{
		if (typeMap.containsKey(gene))
		{
			int s0 = Integer.parseInt(site.substring(1));

			Integer effect = null;
			int closestDist = Integer.MAX_VALUE;

			for (String ss : typeMap.get(gene).keySet())
			{
				Integer eff = typeMap.get(gene).get(ss);

				int s1 = Integer.parseInt(ss.substring(1));

				int dist = Math.abs(s0 - s1);
				if (dist == closestDist && !eff.equals(effect))
				{
					effect = 0;
				}
				else if (dist < closestDist)
				{
					effect = eff;
					closestDist = dist;
				}
			}
			return effect;
		}
		else return null;
	}

	private static List<String> sortSites(Set<String> sites)
	{
		List<String> list = new ArrayList<String>(sites);
		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				try
				{
					return new Integer(o1.substring(1)).compareTo(new Integer(o2.substring(1)));
				}
				catch (NumberFormatException e)
				{
					return 0;
				}
			}
		});
		return list;
	}

	private static List<String> getGenesWithMostSites()
	{
		List<String> genes = new ArrayList<String>(typeMap.keySet());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Integer(typeMap.get(o2).size()).compareTo(typeMap.get(o1).size());
			}
		});
		return genes;
	}

	private static void printSites(String gene)
	{
		System.out.println("Gene: " + gene);
		if (typeMap.containsKey(gene))
		{
			for (String site : sortSites(typeMap.get(gene).keySet()))
			{
				Integer sign = typeMap.get(gene).get(site);
				System.out.print("\tsite: " + site + "\t" + (sign == 1 ? "activating" : sign == -1 ? "inhibiting" : "complex"));
				System.out.println("\t(" + actualMap.get(gene).get(site) + ")");
			}
		}
	}

	static
	{
		typeMap = new HashMap<String, Map<String, Integer>>();
		actualMap = new HashMap<String, Map<String, String>>();
		Scanner sc = new Scanner(PhosphoSitePlus.class.getResourceAsStream("Regulatory_sites"));

		for (int i = 0; i < 4; i++) sc.nextLine();

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] split = line.split("\t");
			if (split.length < 13) continue;

			if (!split[6].equals("human")) continue;

			String gene = HGNC.getSymbol(split[4]);
			if (gene == null) continue;

			if (!split[8].equals("PHOSPHORYLATION")) continue;

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<String, Integer>());
			if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<String, String>());

			String site = split[7];
			actualMap.get(gene).put(site, split[12]);

			boolean actWord = false;
			boolean inhWord = false;

			if ((split[12].contains("induced") &&
				!split[12].contains("receptor desensitization, induced")))
			{
				actWord = true;
			}
			if ((split[12].contains("inhibited") ||
				split[12].contains("receptor desensitization, induced")))
			{
				inhWord = true;
			}

			if (actWord == inhWord)
			{
				if (split[12].contains("stabilization"))
				{
					actWord = true;
				}
				if (split[12].contains("degradation"))
				{
					inhWord = true;
				}
			}

			if (actWord == inhWord)
			{
				typeMap.get(gene).put(site, 0);
			}
			else
			{
				typeMap.get(gene).put(site, actWord ? 1 : -1);
			}
		}

		sc = new Scanner(PhosphoSitePlus.class.getResourceAsStream("manually-curated-sites.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");
			if (token.length < 3) continue;
			String gene = token[0];

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<String, Integer>());
			if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<String, String>());

			String site = token[1];
			int sign = Integer.parseInt(token[2]);

			typeMap.get(gene).put(site, sign);
			actualMap.get(gene).put(site, "manual curation");
		}
	}

	static void printUniqueAA()
	{
		Set<String> sites = new HashSet<String>();
		for (String gene : typeMap.keySet())
		{
			sites.addAll(typeMap.get(gene).keySet());
		}
		TermCounter tc = new TermCounter();
		for (String site : sites)
		{
			tc.addTerm(site.substring(0, 1));
		}
		tc.print();
	}

	public static void main(String[] args)
	{
//		List<String> list = getGenesWithMostSites();
//		for (int i = 0; i < 10; i++)
//		{
//			printSites(list.get(i));
//		}
		printSites("ERBB3");
//		printUniqueAA();
	}
}