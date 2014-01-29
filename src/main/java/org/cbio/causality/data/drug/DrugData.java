package org.cbio.causality.data.drug;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class DrugData
{
	private static Map<String, Set<String>> drug2target = new HashMap<String, Set<String>>();
	private static Map<String, Set<String>> target2drug = new HashMap<String, Set<String>>();
	private static Map<String, String> description = new HashMap<String, String>();

	private static Set<String> fdaApproved = new HashSet<String>();
	private static Set<String> cancerDrug = new HashSet<String>();
	private static Set<String> nutraceutical = new HashSet<String>();

	public static Set<String> getDrugs(String target)
	{
		if (!target2drug.containsKey(target)) return Collections.emptySet();
		return target2drug.get(target);
	}

	public static Set<String> getTargets(String drug)
	{
		if (!drug2target.containsKey(drug)) return Collections.emptySet();
		return drug2target.get(drug);
	}

	public static boolean isFDAApproved(String drug)
	{
		return fdaApproved.contains(drug);
	}

	public static boolean isCancerDrug(String drug)
	{
		return cancerDrug.contains(drug);
	}

	public static boolean isNutraceutical(String drug)
	{
		return nutraceutical.contains(drug);
	}

	public static String getDescription(String drug)
	{
		return description.get(drug);
	}

	public static Set<String> getFDAApprovedDrugs(String target)
	{
		Set<String> drugs = new HashSet<String>(getDrugs(target));
		drugs.retainAll(fdaApproved);
		return drugs;
	}
	
	public static Set<String> getCancerDrugs(String target)
	{
		Set<String> drugs = new HashSet<String>(getDrugs(target));
		drugs.retainAll(cancerDrug);
		return drugs;
	}

	public static Map<String, Set<String>> getDrugs(Set<String> targets)
	{
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();
		for (String target : targets)
		{
			for (String drug : getDrugs(target))
			{
				if (!map.containsKey(drug)) map.put(drug, new HashSet<String>());
				map.get(drug).add(target);
 			}
		}
		return map;
	}

	public static List<String> sortDrugs(final Map<String, Set<String>> drugs)
	{
		List<String> sorted = new ArrayList<String>(drugs.keySet());

		Collections.sort(sorted, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				Integer size1 = drugs.get(o1).size();
				Integer size2 = drugs.get(o2).size();

				if (size1.equals(size2)) return o1.compareTo(o2);
				return size2.compareTo(size1);
			}
		});

		return sorted;
	}

	private static String getHeader(Scanner scan)
	{
		while(scan.hasNextLine())
		{
			String line = scan.nextLine();
			if (!line.startsWith("#")) return line;
		}
		return null;
	}

	static
	{
		Scanner scan = new Scanner(DrugData.class.getResourceAsStream("drugs.tsv"));
		String header = getHeader(scan);
		assert header.equals("PiHelper_Drug_ID\tDrug_Name\tDrug_Synonyms\tDescription\tNumber_of_Targets\tATC_Codes\tisFdaApproved\tisCancerDrug\tisNutraceutical\tNumber_Of_Clinical_Trials\tDataSources\tReferences") : "Something changed in matrix!";

		while (scan.hasNextLine())
		{
			String line = scan.nextLine();
			String[] tok = line.split("\t");

			if (tok.length > 8)
			{
				if (tok[3].length() > 2) description.put(tok[1], tok[3].substring(1, tok[3].length() - 1));
				if (tok[6].startsWith("t")) fdaApproved.add(tok[1]);
				if (tok[7].startsWith("t")) cancerDrug.add(tok[1]);
				if (tok[8].startsWith("t")) nutraceutical.add(tok[1]);
			}
		}

		scan = new Scanner(DrugData.class.getResourceAsStream("drugtargets.tsv"));
		header = getHeader(scan);
		assert header.equals("PiHelper_DrugTarget_ID\tHGNC_Symbol\tDrug_Name\tDataSources\tReferences") : "Something changed in matrix!";

		while (scan.hasNextLine())
		{
			String line = scan.nextLine();
			String[] tok = line.split("\t");

			if (tok.length > 2)
			{
				if (!drug2target.containsKey(tok[2])) drug2target.put(tok[2], new HashSet<String>());
				if (!target2drug.containsKey(tok[1])) target2drug.put(tok[1], new HashSet<String>());

				drug2target.get(tok[2]).add(tok[1]);
				target2drug.get(tok[1]).add(tok[2]);
			}
		}
	}

	public static void main(String[] args)
	{
		String drug = "Cetuximab";
		System.out.println(getTargets(drug));
		System.out.println(isCancerDrug(drug));
		System.out.println(isFDAApproved(drug));
		System.out.println(isNutraceutical(drug));
		System.out.println(getDescription(drug));

		System.out.println(getDrugs("EGFR"));
	}
}
