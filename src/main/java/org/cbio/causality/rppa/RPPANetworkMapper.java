package org.cbio.causality.rppa;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.PhosphoGraph;
import org.cbio.causality.analysis.Relation;
import org.cbio.causality.network.SignedPC;
import org.cbio.causality.signednetwork.SignedType;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class RPPANetworkMapper
{
	private static Map<SignedType,Graph> graph;

	public static List<Relation> map(Collection<RPPAData> data)
	{
		List<Relation> rels = new ArrayList<Relation>();
		if (graph == null) graph = SignedPC.getAllGraphs();

		for (RPPAData d1 : data)
		{
			findRelations(d1, data, rels);
		}
		return rels;
	}

	private static void findRelations(RPPAData d1, Collection<RPPAData> data, List<Relation> rels)
	{
		for (RPPAData d2 : data)
		{
			if (d1 == d2 || d1.genes.equals(d2.genes)) continue;

			for (String source : d1.genes)
			{
				for (String target : d2.genes)
				{
					if (source.equals(target)) continue;

					for (SignedType type : graph.keySet())
					{
						if (graph.get(type).getDownstream(source).contains(target))
						{
							boolean match = false;
							int sign = 0;

							switch (type)
							{
								case PHOSPHORYLATES:
								{
									if (d2.isPhospho())
									{
										match = true;
										sign = d1.getSelfEffect();
									}
								} break;
								case DEPHOSPHORYLATES:
								{
									if (d2.isPhospho())
									{
										match = true;
										sign = -d1.getSelfEffect();
									}
								} break;
								case UPREGULATES_EXPRESSION:
								{
									if (d2.isTotalProt())
									{
										match = true;
										sign = d1.getSelfEffect();
									}
								} break;
								case DOWNREGULATES_EXPRESSION:
								{
									if (d2.isTotalProt())
									{
										match = true;
										sign = -d1.getSelfEffect();
									}
								} break;
							}

							if (match)
							{
								Set<String> sites = (graph.get(type) instanceof PhosphoGraph) ?
									((PhosphoGraph) graph.get(type)).getSites(source, target) :
									null;

								rels.add(new Relation(source, target, type, d1, d2, sign,
									graph.get(type).getMediatorsInString(source, target),
									sites));
							}
						}
					}
				}
			}
		}
	}

	public static void removeConflictingAndInsignificant(List<Relation> rels)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (rel.dataChangesAsUnxpected() || rel.dataChangesInsignificant()) iter.remove();
		}
	}

	public static void keepMatching(List<Relation> rels, boolean siteMatch)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (!rel.dataChangesAsExpected()) iter.remove();
			else if (siteMatch && !rel.siteMatches()) iter.remove();
		}
	}

	public static void keepConflicting(List<Relation> rels, boolean siteMatch)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (!rel.dataChangesAsUnxpected()) iter.remove();
			else if (siteMatch && !rel.siteMatches()) iter.remove();
		}
	}

	public static void keepChangedEndsOnly(List<Relation> rels)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (rel.sourceData.getChangeSign() == 0 || rel.targetData.getChangeSign() == 0)
				iter.remove();
		}
	}

	public static void cropTo(List<Relation> rels, Collection<String> cropTo)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (!cropTo.contains(rel.source) && !cropTo.contains(rel.target)) iter.remove();
		}
	}

	public static void removeConflictingActivities(List<Relation> rels, Collection<RPPAData> datas)
	{
		Map<String, Integer> map = readActivities(datas);

		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (map.containsKey(rel.source) && rel.sourceData.getActvityChangeSign() != map.get(rel.source))
				iter.remove();
		}
	}

	private static final Color MIN_UP = new Color(255, 210, 210);
	private static final Color MIN_DOWN = new Color(210, 210, 255);
	private static final Color MAX_UP = new Color(255, 80, 40);
	private static final Color MAX_DOWN = new Color(40, 80, 255);

	private static final DecimalFormat FMT = new DecimalFormat("0.##");

	public static void writeGraph(Collection<RPPAData> datas, double thresholdVal, String filename,
		GraphType type, Collection<String> cropTo) throws IOException
	{
		checkForDuplicateTotalProt(datas);
		List<Relation> relations = RPPANetworkMapper.map(datas);

		if (cropTo != null && !cropTo.isEmpty()) cropTo(relations, cropTo);

		filterToDesiredSubset(relations, datas, type);

		System.out.println("causative relations = " + relations.size());

		Set<String> nodesInRels = new HashSet<String>();
		Set<String> activated = new HashSet<String>();
		Set<String> inactivated = new HashSet<String>();
		Set<String> sifLines = new HashSet<String>();
		for (Relation rel : relations)
		{
			sifLines.add(rel.getEdgeData());
			nodesInRels.addAll(rel.sourceData.genes);
			nodesInRels.addAll(rel.targetData.genes);
			int sign = rel.sourceData.getActvityChangeSign();
			if (sign == 1) activated.add(rel.source);
			else if (sign == -1) inactivated.add(rel.source);
		}

		BufferedWriter writer;

		if (type != GraphType.EXISTING_NETWORK)
		{
			writer = new BufferedWriter(new FileWriter(filename));
			for (String line : sifLines) writer.write(line + "\n");

			for (RPPAData data : datas)
			{
				if (data.getChangeSign() != 0)
				{
					for (String gene : data.genes)
					{
//						if (!nodesInRels.contains(gene))
						{
							writer.write(gene + "\n");
							nodesInRels.add(gene);
						}
					}
				}
			}

			writer.close();
		}

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		writer = new BufferedWriter(new FileWriter(filename));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");

		for (Relation rel : relations)
		{
			writer.write("edge\t" + rel.source + " " + rel.edgeType.getTag() + " " + rel.target +
				"\tcolor\t" + getEdgeColor(rel.edgeType) + "\n");
		}

		double maxVal = 0;
		for (RPPAData data : datas)
		{
			double val = Math.abs(data.getChangeValue());
			if (val > maxVal) maxVal = val;
		}

		Map<String, RPPAData> totalProtMap = getTotalProtMapping(datas);

		for (RPPAData data : datas)
		{
			for (String gene : data.genes)
			{
				int changeSign = data.getChangeSign();
				if (data.isPhospho() || (data.isTotalProt() && totalProtMap.get(gene) != data))
				{
					String colS = getColor(data.getChangeValue(), thresholdVal, maxVal,
						changeSign > 0 ? MIN_UP : MIN_DOWN, changeSign > 0 ? MAX_UP : MAX_DOWN);
					String bor = data.effect == null || data.effect == RPPAData.SiteEffect.COMPLEX ?
						"0 0 0" : data.effect == RPPAData.SiteEffect.ACTIVATING ? "0 180 20" : "180 0 20";
					String let = data.isPhospho() ? "p" : "t";

					writer.write("node\t" + gene + "\trppasite\t" + data.id + "(" +
						FMT.format(data.getChangeValue()) + ")|" + let + "|" +
						colS + "|" + bor + "\n");
				}
				else if (data.isActivity())
				{
					String colS = "255 255 255";
					String bor = changeSign > 0 ? "0 180 20" : "180 0 20";

					writer.write("node\t" + gene + "\trppasite\t" + data.id + "|" +
						(changeSign > 0 ? "a" : "i") + "|" + colS + "|" + bor + "\n");
				}
				else
				{
					writer.write("node\t" + gene + "\tcolor\t" +
						getColor(data.getChangeValue(), thresholdVal, maxVal,
							changeSign > 0 ? MIN_UP : MIN_DOWN,
							changeSign > 0 ? MAX_UP : MAX_DOWN) + "\n");

					writer.write("node\t" + gene + "\ttooltip\t" + data.id +  "(" +
						FMT.format(data.getChangeValue()) + ")\n");
				}
			}
		}

		for (String name : activated)
		{
			if (inactivated.contains(name))
			{
				writer.write("node\t" + name + "\tbordercolor\t180 100 100\n");
			}
		}

		Map<String, Set<String>> gene2rppa = new HashMap<String, Set<String>>();
		for (RPPAData data : datas)
		{
			for (String gene : data.genes)
			{
				if (!gene2rppa.containsKey(gene)) gene2rppa.put(gene, new HashSet<String>());
				gene2rppa.get(gene).add(data.id);
			}
		}

		writer.close();
	}

	private static void filterToDesiredSubset(List<Relation> relations, Collection<RPPAData> datas, GraphType type)
	{
		switch (type)
		{
			case COMPATIBLE: RPPANetworkMapper.keepMatching(relations, false); break;
			case COMPATIBLE_WITH_SITE_MATCH: RPPANetworkMapper.keepMatching(relations, true); break;
			case CONFLICTING: RPPANetworkMapper.keepConflicting(relations, false); break;
			case CONFLICTING_WITH_SITE_MATCH: RPPANetworkMapper.keepConflicting(relations, true); break;
			case NON_CONFLICTING: RPPANetworkMapper.removeConflictingAndInsignificant(relations); break;
			case CHANGED_ONLY: RPPANetworkMapper.keepChangedEndsOnly(relations); break;
		}

		if (type == GraphType.COMPATIBLE || type == GraphType.COMPATIBLE_WITH_SITE_MATCH)
		{
			removeConflictingActivities(relations, datas);

			// Decides an activity if it is ambiguous, based on majority of compatible relations
//			Set<Relation> underSupported = getUndersupportedRelations(relations,
//				getCropped(datas, findChangingGenesInConflict(getGeneToRPPAMap(datas))));
//			relations.removeAll(underSupported);
		}
	}

	private static String getColor(double val, double minVal, double maxVal, Color minCol, Color maxCol)
	{
		val = Math.abs(val);
		if (val < minVal)
		{
			return "250 250 250";
		}

		double ratio = (val - minVal) / (maxVal - minVal);
		if (ratio > 1) ratio = 1;

		return (
			(int) Math.round(minCol.getRed() + (ratio * (maxCol.getRed() - minCol.getRed()))) + " " +
				(int) Math.round(minCol.getGreen() + (ratio * (maxCol.getGreen() - minCol.getGreen()))) + " " +
				(int) Math.round(minCol.getBlue() + (ratio * (maxCol.getBlue() - minCol.getBlue()))));
	}

	private static String getEdgeColor(SignedType type)
	{
		switch (type)
		{
			case PHOSPHORYLATES:
			case UPREGULATES_EXPRESSION: return "0 100 0";
			case DEPHOSPHORYLATES:
			case DOWNREGULATES_EXPRESSION: return "100 0 0";
			default: return null;
		}
	}

	// Section: selecting most likely activity by looking at the changes at the downstream.

	protected static Map<String, Set<RPPAData>> getGeneToRPPAMap(Collection<RPPAData> datas)
	{
		Map<String, Set<RPPAData>> map = new HashMap<String, Set<RPPAData>>();
		for (RPPAData data : datas)
		{
			for (String gene : data.genes)
			{
				if (!map.containsKey(gene)) map.put(gene, new HashSet<RPPAData>());
				map.get(gene).add(data);
			}
		}
		return map;
	}

	protected static Set<String> findChangingGenesInConflict(Map<String, Set<RPPAData>> map)
	{
		Set<String> set = new HashSet<String>();

		for (String gene : map.keySet())
		{
			boolean up = false;
			boolean down = false;

			for (RPPAData data : map.get(gene))
			{
				int sign = data.getActvityChangeSign();
				if (sign > 0) up = true;
				else if (sign < 0) down = true;
			}

			if (up && down) set.add(gene);
		}
		return set;
	}

	protected static Set<RPPAData> getCropped(Collection<RPPAData> datas, Set<String> genes)
	{
		Set<RPPAData> set = new HashSet<RPPAData>();
		for (RPPAData data : datas)
		{
			for (String gene : data.genes)
			{
				if (genes.contains(gene))
				{
					set.add(data);
					break;
				}
			}
		}
		return set;
	}

	protected static Map<String, Integer> readActivities(Collection<RPPAData> datas)
	{
		Map<String, Integer> map = new HashMap<String, Integer>();

		for (RPPAData data : datas)
		{
			if (data.isActivity())
			{
				for (String gene : data.genes)
				{
					map.put(gene, data.getActvityChangeSign());
				}
			}
		}
		return map;
	}

	protected static Set<Relation> getUndersupportedRelations(Collection<Relation> rels, Set<RPPAData> datas)
	{
		Map<String, Set<String>> actSupport = new HashMap<String, Set<String>>();
		Map<String, Set<String>> inaSupport = new HashMap<String, Set<String>>();

		Map<String, Set<Relation>> actRels = new HashMap<String, Set<Relation>>();
		Map<String, Set<Relation>> inaRels = new HashMap<String, Set<Relation>>();

		for (Relation rel : rels)
		{
			if (!datas.contains(rel.sourceData)) continue;

			Map<String, Set<String>> map = null;
			Map<String, Set<Relation>> relMap = null;

			if (rel.dataChangesAsExpected() && rel.siteMatches())
			{
				if (rel.sourceData.getActvityChangeSign() > 0)
				{
					map = actSupport;
					relMap = actRels;
				}
				else if (rel.sourceData.getActvityChangeSign() < 0)
				{
					map = inaSupport;
					relMap = inaRels;
				}
			}

			if (map != null)
			{
				if (!map.containsKey(rel.source)) map.put(rel.source, new HashSet<String>());
				map.get(rel.source).add(rel.target);

				if (!relMap.containsKey(rel.source)) relMap.put(rel.source, new HashSet<Relation>());
				relMap.get(rel.source).add(rel);
			}
		}

		Set<String> genes = new HashSet<String>(actSupport.keySet());
		genes.addAll(inaSupport.keySet());

		Set<Relation> under = new HashSet<Relation>();

		for (String gene : genes)
		{
			int actCnt = actSupport.containsKey(gene) ? actSupport.get(gene).size() : 0;
			int inaCnt = inaSupport.containsKey(gene) ? inaSupport.get(gene).size() : 0;
			System.out.println("gene = " + gene + "\tact = " + actCnt + "\tina = " + inaCnt);

			if (actCnt > 0 && actCnt < inaCnt) under.addAll(actRels.get(gene));
			else if (inaCnt > 0 && actCnt > inaCnt) under.addAll(inaRels.get(gene));
		}

		return under;
	}

	protected static void checkForDuplicateTotalProt(Collection<RPPAData> datas)
	{
		Map<String, Set<RPPAData>> map = new HashMap<String, Set<RPPAData>>();
		for (RPPAData data : datas)
		{
			if (!data.isTotalProt()) continue;
			for (String gene : data.genes)
			{
				if (!map.containsKey(gene)) map.put(gene, new HashSet<RPPAData>());
				map.get(gene).add(data);
			}
		}
		for (String gene : map.keySet())
		{
			Set<RPPAData> set = map.get(gene);
			if (set.size() > 1)
			{
				System.out.print("Multiple total protein for sym " + gene + ":");
				for (RPPAData d : set)
				{
					System.out.print("    " + d.id);
				}
				System.out.println();
			}
		}
	}

	protected static Map<String, RPPAData> getTotalProtMapping(Collection<RPPAData> datas)
	{
		Map<String, RPPAData> map = new HashMap<String, RPPAData>();
		for (RPPAData data : datas)
		{
			if (!data.isTotalProt()) continue;
			for (String gene : data.genes)
			{
				if (!map.containsKey(gene)) map.put(gene, data);
			}
		}
		return map;
	}

	// Activity / inactivity support

	public static Map<String, int[]> getDownstreamActivitySupportCounts(Collection<RPPAData> data,
		GraphType type)
	{
		if (graph == null) graph = SignedPC.getAllGraphs();

		Map<String, int[]> count = new HashMap<String, int[]>();
		Set<String> symbols = new HashSet<String>();
		for (RPPAData d : data)
		{
			symbols.addAll(d.genes);
		}

		List<Relation> relations = new ArrayList<Relation>();
		for (String symbol : symbols)
		{
			relations.clear();
			RPPAData d = new RPPAData(symbol + "-act", null, Arrays.asList(symbol), null);
			d.makeActivityNode(true);
			findRelations(d, data, relations);
			filterToDesiredSubset(relations, data, type);
			int act = relations.size();

			relations.clear();
			d = new RPPAData(symbol + "-inh", null, Arrays.asList(symbol), null);
			d.makeActivityNode(false);
			findRelations(d, data, relations);
			filterToDesiredSubset(relations, data, type);
			int inh = relations.size();

			count.put(symbol, new int[]{act, inh});
		}

		for (String symbol : count.keySet())
		{
			System.out.println(symbol + "\t" + count.get(symbol)[0] + "\t" + count.get(symbol)[1]);
		}
		return count;
	}


	// Graph size statistics

	public static List<Integer> getNullGraphSizes(List<RPPAData> datas, int iteration,
		GraphType type)
	{
		List<Integer> sizes = new ArrayList<Integer>();

		for (int i = 0; i < iteration; i++)
		{
			RPPAData.shuffleValues(datas);
			List<Relation> relations = map(datas);

			filterToDesiredSubset(relations, datas, type);

			sizes.add(relations.size());
		}
		Collections.sort(sizes);
		Collections.reverse(sizes);
		return sizes;
	}

	public static enum GraphType
	{
		EXISTING_NETWORK,
		ALL_INCLUSIVE,
		CHANGED_ONLY,
		NON_CONFLICTING,
		COMPATIBLE,
		COMPATIBLE_WITH_SITE_MATCH,
		CONFLICTING,
		CONFLICTING_WITH_SITE_MATCH
	}
}
