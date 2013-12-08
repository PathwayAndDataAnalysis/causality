package org.cbio.causality.ubique;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Searcher;
import org.biopax.paxtools.pattern.miner.ChemicalAffectsThroughBindingMiner;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.biopax.paxtools.pattern.util.HGNC;
import org.cbio.causality.analysis.Traverse;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class UbiqueDetection
{
	SimpleIOHandler handler = new SimpleIOHandler();

	public static void main(String[] args) throws IOException
	{
		UbiqueDetection det = new UbiqueDetection();
		det.printStats();
	}

	private void printStats() throws FileNotFoundException
	{
		Traverse trav = new Traverse();
		trav.load(SIFProducer.SIF_FILENAME, Collections.<String>emptySet(), new HashSet<String>(
			Arrays.asList(SIFType.CONSUMPTION_CONTROLLED_BY.getTag(),
				SIFType.CONTROLS_PRODUCTION_OF.getTag())));

		// calculate degrees

		final Map<String, Integer> degreeMap = new HashMap<String, Integer>();
		Map<String, Integer> indegreeMap = new HashMap<String, Integer>();
		Map<String, Integer> outdegreeMap = new HashMap<String, Integer>();

		Map<String, Set<String>> upstrMap = new HashMap<String, Set<String>>();
		Map<String, Set<String>> dwstrMap = new HashMap<String, Set<String>>();
		Map<String, Set<String>> neighMap = new HashMap<String, Set<String>>();

		for (String name : trav.getSymbols())
		{
			if (HGNC.getSymbol(name) != null) continue;

			Set<String> upstr = trav.goBFS(name, false);
			Set<String> dwstr = trav.goBFS(name, true);

			neighMap.put(name, new HashSet<String>());
			neighMap.get(name).addAll(upstr);
			neighMap.get(name).addAll(dwstr);

//			Set<String> temp = new HashSet<String>(upstr);
//			upstr.removeAll(dwstr);
//			dwstr.removeAll(temp);

			int indegree = upstr.size();
			int outdegree = dwstr.size();

			int degree = indegree + outdegree;

			if (degree > 0)
			{
				degreeMap.put(name, neighMap.get(name).size());
				indegreeMap.put(name, indegree);
				outdegreeMap.put(name, outdegree);

				upstrMap.put(name, upstr);
				dwstrMap.put(name, dwstr);
			}
		}

		// calculate clustering coefficients

		trav.clear();
		trav.load(SIFProducer.SIF_FILENAME, Collections.singleton(SIFType.NEIGHBOR_OF.getTag()),
				Collections.<String>emptySet());

		Map<String, Double> inCluster = new HashMap<String, Double>();
		Map<String, Double> outCluster = new HashMap<String, Double>();
		Map<String, Double> cluster = new HashMap<String, Double>();

		for (String name : degreeMap.keySet())
		{
			inCluster.put(name, calcClusterCoef(upstrMap.get(name), trav));
			outCluster.put(name, calcClusterCoef(dwstrMap.get(name), trav));
			cluster.put(name, calcClusterCoef(neighMap.get(name), trav));
		}


		Model model = handler.convertFromOWL(new FileInputStream(
				"/home/ozgun/Projects/biopax-pattern/All-Human-Data.owl"));

		Set<BiochemicalReaction> inters = collectInteractionsWithChemicalParticipantOnly(model);
		System.out.println("interaction size = " + inters.size());
		Map<String, Integer> buddyCounts = getBuddyCounts(inters);
		Map<String, Integer> minorCounts = getUbiqueScores(inters, buddyCounts);

		Map<String, Integer> interCnt = getParticipantOfCnt(model, minorCounts.keySet(), inters);
		Map<String, Integer> controlCnt = getControllerOfCnt(model, minorCounts.keySet());

		Map<String, Integer> leftEssentCnt = getEssentialInterCounts(inters, buddyCounts, true);
		Map<String, Integer> rightEssentCnt = getEssentialInterCounts(inters, buddyCounts, false);
		Map<String, Integer> bothSideCnt = getBothSideCounts(inters);

		// print

		List<String> names = new ArrayList<String>(degreeMap.keySet());
		Collections.sort(names, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return degreeMap.get(o2).compareTo(degreeMap.get(o1));
			}
		});

		System.out.print("Name\tdegree\toverall clustering\tconsumed\tconsumer clustering\t" +
				"produced\tproducer clustering\tminors\tbuddies\tedges\tcontrols\tleft " +
				"essen\tright essen\tboth-side");
		for (String name : names)
		{
			System.out.print("\n" + name + "\t" + degreeMap.get(name) + "\t" +
					cluster.get(name) + "\t" +
					outdegreeMap.get(name) + "\t" + outCluster.get(name) + "\t" +
					indegreeMap.get(name) + "\t" + inCluster.get(name) + "\t" +
					minorCounts.get(name) + "\t" + buddyCounts.get(name) + "\t" +
					interCnt.get(name) + "\t" + controlCnt.get(name) + "\t" +
					leftEssentCnt.get(name) + "\t" + rightEssentCnt.get(name) + "\t" +
					bothSideCnt.get(name)
			);
		}
	}

	private double calcClusterCoef(Set<String> genes, Traverse trav)
	{
		if (genes.isEmpty()) return Double.NaN;
		else if (genes.size() == 1) return 1;

		int foundEdges = 0;
		for (String gene : genes)
		{
			Set<String> neighbors = trav.getNeighbors(gene);
			neighbors.remove(gene);
			neighbors.retainAll(genes);
			foundEdges += neighbors.size();
		}

		double coef = foundEdges / (double) (genes.size() * (genes.size() - 1));

		return coef;
	}

	public void scoreUbiques() throws IOException
	{
		Model model = handler.convertFromOWL(new FileInputStream(
			"/home/ozgun/Projects/biopax-pattern/All-Human-Data.owl"));

		Set<BiochemicalReaction> inters = collectInteractionsWithChemicalParticipantOnly(model);
		System.out.println("interaction size = " + inters.size());
		Map<String, Integer> buddyCounts = getBuddyCounts(inters);
		Map<String, Integer> minorCounts = getUbiqueScores(inters, buddyCounts);

		Map<String, Integer> interCnt = getParticipantOfCnt(model, minorCounts.keySet(), inters);
		Map<String, Integer> controlCnt = getControllerOfCnt(model, minorCounts.keySet());

		Map<String, Integer> leftEssentCnt = getEssentialInterCounts(inters, buddyCounts, true);
		Map<String, Integer> rightEssentCnt = getEssentialInterCounts(inters, buddyCounts, false);
		Map<String, Integer> bothSideCnt = getBothSideCounts(inters);

		Map<String, Integer> geneCounts = getGeneCounts(model, minorCounts.keySet());

//		Map<String, Set<String>> minors = getMinors(inters, buddyCounts);

		model = null;

		Map<String, Double> clusCoef = getClusteringCoefficients(minorCounts.keySet());

		List<String> sorted = getSortedList(minorCounts);

		System.out.print("Name\tk.safdas\tsafdas\tedges\tctrls\tleft essent.\tright essent.\tbothside\tbinding prot\tcluster coef");

		for (String name : sorted)
		{
			System.out.print("\n" + name + "\t" + minorCounts.get(name) + "\t" +
				buddyCounts.get(name) + "\t" + interCnt.get(name) + "\t" + controlCnt.get(name) +
				"\t" + (leftEssentCnt.containsKey(name) ? leftEssentCnt.get(name) : 0) + "\t" +
				(rightEssentCnt.containsKey(name) ? rightEssentCnt.get(name) : 0) + "\t" +
				(bothSideCnt.containsKey(name) ? bothSideCnt.get(name) : 0) + "\t" +
				(geneCounts.containsKey(name) ? geneCounts.get(name) : 0) + "\t" +
				clusCoef.get(name));

//			int i = 0;
//			for (String s : minors.get(name))
//			{
//				System.out.print("\t" + s);
//				if (i++ > 10) break;
//			}
		}
	}

	private Set<BiochemicalReaction> collectInteractionsWithChemicalParticipantOnly(Model model)
	{
		Set<BiochemicalReaction> set = new HashSet<BiochemicalReaction>();
		for (BiochemicalReaction inter : model.getObjects(BiochemicalReaction.class))
		{
			int cnt  = 0;
			for (Entity ent : inter.getParticipant())
			{
				if (ent instanceof SmallMolecule) cnt++;
				else
				{
					cnt = 0;
					break;
				}
			}
			if (cnt >= 1)
			{
				set.add(inter);
			}
		}
		return set;
	}

	private Map<String, Integer> getBuddyCounts(Set<BiochemicalReaction> inters)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();
		Map<String, Set<String>> sets = new HashMap<String, Set<String>>();

		for (BiochemicalReaction inter : inters)
		{
			countBuddies(sets, getERs(inter, true));
			countBuddies(sets, getERs(inter, false));
		}

		for (String name : sets.keySet())
		{
			cnt.put(name, sets.get(name).size());
		}

		return cnt;
	}

	private void countBuddies(Map<String, Set<String>> sets, Set<EntityReference> ers)
	{
		Set<String> names = getNames(ers);
		for (String name : names)
		{
			if (!sets.containsKey(name)) sets.put(name, new HashSet<String>());

			sets.get(name).addAll(names);
			sets.get(name).remove(name);
		}
	}

	private Map<String, Integer> getUbiqueScores(Set<BiochemicalReaction> inters,
		Map<String, Integer> cnt)
	{
		Map<String, Integer> score = new HashMap<String, Integer>();
		Map<String, Set<String>> minors = getMinors(inters, cnt);

		for (String name : minors.keySet())
		{
			score.put(name, minors.get(name).size());
		}

		return score;
	}

	private Map<String, Set<String>> getMinors(Set<BiochemicalReaction> inters, Map<String, Integer> cnt)
	{
		Map<String, Set<String>> minors = new HashMap<String, Set<String>>();

		for (BiochemicalReaction inter : inters)
		{
			countMinors(cnt, minors, getERs(inter, true));
			countMinors(cnt, minors, getERs(inter, false));
		}
		return minors;
	}

	private void countMinors(Map<String, Integer> cnt, Map<String, Set<String>> buddies,
		Set<EntityReference> ers)
	{
		if (ers.size() < 2) return;

		Set<String> names = getNames(ers);

		for (String name : names)
		{
			if (!buddies.containsKey(name)) buddies.put(name, new HashSet<String>());

			for (String nm : names)
			{
				if (cnt.get(nm) < cnt.get(name)) buddies.get(name).add(nm);
			}
		}
	}


	private String name(EntityReference er)
	{
		return name(er.getDisplayName());
	}

	private String name(String name)
	{
		return ChemicalNameNormalizer.nornmalize(name);
	}

	private Set<String> getNames(Set<EntityReference> ers)
	{
		Set<String> names = new HashSet<String>();
		for (EntityReference er : ers)
		{
			names.add(name(er));
		}
		return names;
	}

	private Set<EntityReference> getERs(BiochemicalReaction inter, boolean left)
	{
		Set<EntityReference> ers = new HashSet<EntityReference>();
		for (Entity ent : left ? inter.getLeft() : inter.getRight())
		{
			if (ent instanceof SmallMolecule)
			{
				EntityReference er = ((SimplePhysicalEntity) ent).getEntityReference();

				if (er != null) ers.add(er);
			}
		}
		return ers;
	}

	private int getMinCnt(Set<EntityReference> set, Map<String, Integer> cnt)
	{
		int c = Integer.MAX_VALUE;
		for (EntityReference er : set)
		{
			String name = name(er);
			if (cnt.get(name) < c) c = cnt.get(name);
		}
		return c;
	}

	private Set<EntityReference> getERs(Interaction inter)
	{
		Set<EntityReference> ers = new HashSet<EntityReference>();
		for (Entity ent : inter.getParticipant())
		{
			if (ent instanceof SmallMolecule)
			{
				EntityReference er = ((SimplePhysicalEntity) ent).getEntityReference();

				if (er != null) ers.add(er);
			}
		}
		return ers;
	}
	private List<String> getSortedList(final Map<String, Integer> scores)
	{
		List<String> list = new ArrayList<String>(scores.keySet());

		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return scores.get(o2).compareTo(scores.get(o1));
			}
		});

		return list;
	}

	private static final PathAccessor smrToInterAcc = new PathAccessor("SmallMoleculeReference/entityReferenceOf/participantOf:Conversion");
	private static final PathAccessor smrToControlAcc = new PathAccessor("SmallMoleculeReference/entityReferenceOf/controllerOf");

	private Map<String, Integer> getParticipantOfCnt(Model model, Set<String> names,
		Set<BiochemicalReaction> inters)
	{
		return getCounts(model, names, smrToInterAcc, inters);
	}

	private Map<String, Integer> getControllerOfCnt(Model model, Set<String> names)
	{
		return getCounts(model, names, smrToControlAcc, null);
	}

	private Map<String, Integer> getCounts(Model model, Set<String> names, PathAccessor pa,
		Set<BiochemicalReaction> inters)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();

		Map<String, Set<Object>> items = new HashMap<String, Set<Object>>();

		for (SmallMoleculeReference smr : model.getObjects(SmallMoleculeReference.class))
		{
			String name = name(smr);
			if (!names.contains(name)) continue;

			if (!items.containsKey(name)) items.put(name, new HashSet<Object>());

			Set vals = pa.getValueFromBean(smr);
			if (inters != null) vals.retainAll(inters);
			items.get(name).addAll(vals);
		}

		for (String name : items.keySet())
		{
			cnt.put(name, items.get(name).size());
		}
		return cnt;
	}

	private Map<String, Integer> getEssentialInterCounts(Set<BiochemicalReaction> inters,
		Map<String, Integer> degrees, boolean left)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();

		for (BiochemicalReaction inter : inters)
		{
			countEssential(degrees, cnt, getERs(inter, left));
		}

		return cnt;
	}

	private void countEssential(Map<String, Integer> degrees, Map<String, Integer> cnt,
		Set<EntityReference> ers)
	{
		int minCnt = getMinCnt(ers, degrees);
		for (EntityReference er : ers)
		{
			String name = name(er);
			if (!cnt.containsKey(name)) cnt.put(name, 0);

			if (degrees.get(name) == minCnt)
			{
				cnt.put(name, cnt.get(name) + 1);
			}
		}
	}


	private Map<String, Integer> getBothSideCounts(Set<BiochemicalReaction> inters)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();

		for (BiochemicalReaction inter : inters)
		{
			Set<EntityReference> lefts = getERs(inter, true);
			Set<EntityReference> rights = getERs(inter, false);

			lefts.retainAll(rights);

			for (EntityReference er : lefts)
			{
				String name = name(er);
				if (!cnt.containsKey(name)) cnt.put(name, 0);
				cnt.put(name, cnt.get(name) + 1);
			}
		}

		return cnt;
	}


	private Map<String, Integer> getGeneCounts(Model model, Set<String> names)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();
		Map<String, Set<ProteinReference>> map = new HashMap<String, Set<ProteinReference>>();

		ChemicalAffectsThroughBindingMiner miner = new ChemicalAffectsThroughBindingMiner(null);
		Map<BioPAXElement, List<Match>> matches = Searcher.search(model, miner.getPattern());

		for (BioPAXElement ele : matches.keySet())
		{
			SmallMoleculeReference smr = (SmallMoleculeReference) ele;
			String name = name(smr);

			if (names.contains(name))
			{
				if (!map.containsKey(name)) map.put(name, new HashSet<ProteinReference>());

				for (Match match : matches.get(ele))
				{
					ProteinReference pr = (ProteinReference) match.getLast();
					map.get(name).add(pr);
				}
			}
		}

		for (String name : map.keySet())
		{
			cnt.put(name, map.get(name).size());
		}
		return cnt;
	}

	private Map<String, Double> getClusteringCoefficients(Set<String> names)
	{
		Traverse trav = new Traverse();
		trav.load("/home/ozgun/Projects/biopax-pattern/neighbor-of-with-prot-sm.txt",
			Collections.singleton(SIFType.NEIGHBOR_OF.getTag()), Collections.<String>emptySet());

		Map<String, Set<String>> net = new HashMap<String, Set<String>>();

		for (String sym : trav.getSymbols())
		{
			String name = name(sym);

			if (names.contains(name))
			{
				if (!net.containsKey(name)) net.put(name, new HashSet<String>());

				for (String sym2 : trav.getNeighbors(sym))
				{
					String name2 = name(sym2);
					if (names.contains(name2)) sym2 = name2;
					net.get(name).add(sym2);
				}
			}
		}

		Set<String> seconds = new HashSet<String>();
		for (Set<String> set : net.values())
		{
			seconds.addAll(set);
		}

		for (String second : seconds)
		{
			String name = name(second);
			if (names.contains(name)) continue;

			net.put(second, new HashSet<String>());

			for (String neigh : trav.getNeighbors(second))
			{
				name = name(neigh);
				if (names.contains(name)) neigh = name;
				net.get(second).add(neigh);
			}
		}

		System.out.println("net.size() = " + net.size());
		Map<String, Double> coefs = new HashMap<String, Double>();

		for (String name : names)
		{
			if (!net.containsKey(name))
			{
				System.out.println(name + " not in network");
				continue;
			}

			int edgeCnt = 0;
			for (String neigh : net.get(name))
			{
				Set<String> inNeigh = new HashSet<String>(net.get(neigh));
				inNeigh.retainAll(net.get(name));
				edgeCnt += inNeigh.size();
			}
			int possibleCnt = net.get(name).size() * (net.get(name).size() - 1);

			coefs.put(name, edgeCnt / (double) possibleCnt);
		}
		return coefs;
	}
}
