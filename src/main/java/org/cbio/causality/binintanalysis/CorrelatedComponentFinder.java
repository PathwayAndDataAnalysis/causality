package org.cbio.causality.binintanalysis;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.ComponentSorter;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class CorrelatedComponentFinder
{
	private Graph graph;
	private Dataset1 dataset;
	CBioPortalAccessor portalAcc;
	ExpDataManager expMan;
	Map<String, Map<String, Double>> correlations;
	Map<String, Map<String, Double>> expSigPvals;
	Map<String, Map<String, Double>> enrichSigPvals;
	Map<String, List<Set<String>>> compsPos;
	Map<String, List<Set<String>>> compsNeg;
	double fdrThr;
	int componentSizeThr;

	private static final int MIN_GROUP_SIZE = 1;

	public static void main(String[] args) throws IOException
	{
		CorrelatedComponentFinder finder = new CorrelatedComponentFinder(
			Dataset1.UCEC, 0.05, 2);

		finder.calcCorrelations();
		finder.calcComponents();
		finder.writeComponents();
	}

	public CorrelatedComponentFinder(Dataset1 dataset, double fdrThr,
		int componentSizeThr)
		throws IOException
	{
		this.dataset = dataset;
		this.fdrThr = fdrThr;
		this.componentSizeThr = componentSizeThr;
		loadNetwork();
		loadData();
	}

	private void loadNetwork()
	{
		graph = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
//		graph = HPRD.getGraph();
	}

	private void loadData() throws IOException
	{
		portalAcc = new CBioPortalAccessor(dataset.mutCnCallExpZ);
		expMan = new ExpDataManager(portalAcc.getGeneticProfileById(dataset.exp.getProfileID()[0]),
			portalAcc.getCaseListById(dataset.exp.getCaseListID()));
	}

	private static final PearsonsCorrelation COR = new PearsonsCorrelation();

	private double calcCorrelation(String g1, String g2)
	{
		double[] e1 = expMan.get(g1);
		double[] e2 = expMan.get(g2);

		return COR.correlation(e1, e2);
	}

	public void calcCorrelations()
	{
		correlations = new HashMap<String, Map<String, Double>>();
		expSigPvals = new HashMap<String, Map<String, Double>>();

		for (String g1 : graph.getSymbols())
		{
			if (!correlations.containsKey(g1)) correlations.put(g1, new HashMap<String, Double>());

			for (String g2 : graph.getNeighbors(g1))
			{
				if (correlations.get(g1).containsKey(g2)) continue;

				double cor = calcCorrelation(g1, g2);

				if (Double.isNaN(cor)) continue;

				correlations.get(g1).put(g2, cor);
				if (!correlations.containsKey(g2)) correlations.put(g2, new HashMap<String, Double>());
				correlations.get(g2).put(g1, cor);
			}
		}
	}

	public void calcComponents()
	{
		compsPos = new HashMap<String, List<Set<String>>>();
		compsNeg = new HashMap<String, List<Set<String>>>();

		for (String mut : enrichSigPvals.keySet())
		{
			Set<String> up = new HashSet<String>();
			Set<String> dw = new HashSet<String>();

			for (String target : enrichSigPvals.get(mut).keySet())
			{
			}

			ComponentSorter sorter = new ComponentSorter(up, graph);
			List<Set<String>> components = sorter.getComponents(componentSizeThr);
			if (!components.isEmpty()) compsPos.put(mut, components);

			sorter = new ComponentSorter(dw, graph);
			components = sorter.getComponents(componentSizeThr);
			if (!components.isEmpty()) compsNeg.put(mut, components);
		}

		Set<String> resultGeneratingMuts = new HashSet<String>(compsPos.keySet());
		resultGeneratingMuts.addAll(compsNeg.keySet());
		System.out.println("Result generated for " + resultGeneratingMuts.size() + " mutators");
	}

	public void writeComponents() throws IOException
	{
		writeComponents(compsPos, "-up");
		writeComponents(compsNeg, "-dw");
	}

	public void writeComponents(Map<String, List<Set<String>>> comps, String add) throws IOException
	{
		String dir = "temp/components/";
		new File(dir).mkdirs();
		for (String mut : comps.keySet())
		{
			if (comps.get(mut).isEmpty()) continue;

			BufferedWriter w1 = new BufferedWriter(new FileWriter(dir + mut + add + ".sif"));
			BufferedWriter w2 = new BufferedWriter(new FileWriter(dir + mut + add + ".format"));

			for (Set<String> comp : comps.get(mut))
			{
				for (String g1 : comp)
				{
					Set<String> neighs = graph.isDirected() ?
						graph.getDownstream(g1) : graph.getNeighbors(g1);

					for (String g2 : comp)
					{
						if (graph.isUndirected() && g2.compareTo(g1) < 0) continue;

						if (neighs.contains(g2))
						{
							w1.write(g1 + "\t" + graph.getEdgeType() + "\t" + g2 + "\n");
						}
					}

//					w2.write("node\t" + g1 + "\tcolor\t" + val2Color(correlations.get(mut).get(g1),
//						targetChange.get(mut).get(g1) ? 1 : -1) + "\n");
				}
			}

			w1.close();
			w2.close();
		}
	}

	private String val2Color(double pval, int type)
	{
		double score = Math.min(10, -Math.log(pval));

		int v = 255 - (int) Math.round((255D / 10D) * score);

		switch (type)
		{
			case 1: return "255 " + v + " " + v;
			case -1: return v + " 255 " + v;
			default: throw new IllegalArgumentException();
		}
	}

}
