package org.cbio.causality.network;

import org.cbio.causality.idmapping.HGNC;

import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class PhosphoSitePlus
{
	static Map<String, Map<String, Boolean>> typeMap;

	public static Boolean getEffect(String gene, String site)
	{
		if (typeMap.containsKey(gene))
		{
			return typeMap.get(gene).get(site);
		}
		return null;
	}

	public static Boolean getMajorityEffect(String gene)
	{
		if (typeMap.containsKey(gene))
		{
			int cp = 0;
			int cn = 0;

			for (String site : typeMap.get(gene).keySet())
			{
				if (typeMap.get(gene).get(site)) cp++; else cn++;
			}

			if (cp == cn) return null;

			return cp > cn;
		}
		return null;
	}

	static
	{
		typeMap = new HashMap<String, Map<String, Boolean>>();
		Scanner sc = new Scanner(HPRD.class.getResourceAsStream("Regulatory_sites"));

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

			if (!split[12].contains("induced") && !split[12].contains("inhibited")) continue;
			if (split[12].contains("induced") && split[12].contains("inhibited")) continue;

			String site = split[7];

			boolean sign = split[12].contains("induced");

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<String, Boolean>());

			typeMap.get(gene).put(site, sign);
		}
	}

	public static void main(String[] args)
	{
		System.out.println("getMajorityEffect = " + getMajorityEffect("RPS6"));
	}
}
