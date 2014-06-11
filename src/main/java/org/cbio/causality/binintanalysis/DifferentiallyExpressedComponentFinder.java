package org.cbio.causality.binintanalysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.ComponentSorter;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.BroadAccessor;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.network.HPRD;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;
import org.cbio.causality.util.StudentsT;
import org.cbio.causality.util.Summary;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class DifferentiallyExpressedComponentFinder
{
	private Graph graph;
	private Dataset1 dataset;
	CBioPortalAccessor portalAcc;
	ExpDataManager expMan;
	Map<String, Map<String, Double>> expPvals;
	Map<String, Map<String, Double>> expSigPvals;
	Map<String, Map<String, Double>> enrichSigPvals;
	Map<String, Map<String, Boolean>> targetChange;
	Map<String, List<Set<String>>> compsUp;
	Map<String, List<Set<String>>> compsDw;
	Set<String> mutsig;
	double mutsigThr;
	double fdrThr;
	int componentSizeThr;

	private static final int MIN_GROUP_SIZE = 1;

	public static void main(String[] args) throws IOException
	{
		DifferentiallyExpressedComponentFinder finder = new DifferentiallyExpressedComponentFinder(
			Dataset1.UCEC, 0.05, 0.05, 2);

		finder.calcDiffExpPvals();
		finder.calcComponents();
		finder.writeComponents();
	}

	public DifferentiallyExpressedComponentFinder(Dataset1 dataset, double mutsigThr, double fdrThr,
		int componentSizeThr)
		throws IOException
	{
		this.dataset = dataset;
		this.mutsigThr = mutsigThr;
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
		expMan = new ExpDataManager(portalAcc.getGeneticProfileById(dataset.exp.profileID[0]),
			portalAcc.getCaseListById(dataset.exp.caseListID));

		mutsig = BroadAccessor.getMutsigGenes(dataset.code(), mutsigThr, true);
		removeMutsigWithMissingData();
		System.out.println("mutsig in network = " + mutsig.size());
	}

	private void removeMutsigWithMissingData()
	{
		Iterator<String> iter = mutsig.iterator();
		while (iter.hasNext())
		{
			String ms = iter.next();

			AlterationPack alts = portalAcc.getAlterations(ms);

			if (alts == null) iter.remove();
			else if (alts.get(Alteration.MUTATION) == null ||
				alts.get(Alteration.COPY_NUMBER) == null ||
				alts.get(Alteration.EXPRESSION) == null) iter.remove();
		}
	}


	private double calcDiffPval(String symbol, boolean[] set1, boolean[] set2)
	{
		double[] exp = expMan.get(symbol);

		if (exp == null) return Double.NaN;

		boolean[] cnUnchanged = getCopyNumberUnchanged(symbol);
		if (cnUnchanged == null) return Double.NaN;

		double[] vals1 = getSubset(exp, set1, cnUnchanged);
		double[] vals2 = getSubset(exp, set2, cnUnchanged);

		if (vals1.length < MIN_GROUP_SIZE || vals2.length < MIN_GROUP_SIZE) return Double.NaN;

		return StudentsT.getPValOfMeanDifference(vals1, vals2);
	}

	private double calcMeanChange(String symbol, boolean[] set1, boolean[] set2)
	{
		double[] exp = expMan.get(symbol);

		if (exp == null) return Double.NaN;

		boolean[] cnUnchanged = getCopyNumberUnchanged(symbol);
		if (cnUnchanged == null) return Double.NaN;

		double[] vals1 = getSubset(exp, set1, cnUnchanged);
		double[] vals2 = getSubset(exp, set2, cnUnchanged);

		if (vals1.length < MIN_GROUP_SIZE || vals2.length < MIN_GROUP_SIZE) return Double.NaN;

		return Summary.calcChangeOfMean(vals1, vals2);
	}

	private double[] getSubset(double[] vals, boolean[] inds, boolean[] cnUnchanged)
	{
		List<Double> list = new ArrayList<Double>(vals.length);

		for (int i = 0; i < vals.length; i++)
		{
			if (inds[i] && cnUnchanged[i] && !Double.isNaN(vals[i])) list.add(vals[i]);
		}
		double[] sub = new double[list.size()];
		for (int i = 0; i < sub.length; i++)
		{
			sub[i] = list.get(i);
		}
		return sub;
	}

	private boolean[] getCopyNumberUnchanged(String symbol)
	{
		AlterationPack alts = portalAcc.getAlterations(symbol);

		if (alts == null) return null;

		Change[] cnc = alts.get(Alteration.COPY_NUMBER);

		if (cnc == null) return null;

		boolean[] b = new boolean[cnc.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !cnc[i].isAltered() && !cnc[i].isAbsent();
		}
		return b;
	}

	/**
	 * Gets the samples that the gene is either mutated or intact. In both cases it is expressed,
	 * and not copy number lost.
	 */
	private boolean[] selectSubset(String symbol, boolean mutated)
	{
		AlterationPack alts = portalAcc.getAlterations(symbol);

		if (alts == null) return null;

		Change[] muts = alts.get(Alteration.MUTATION);
		Change[] cncs = alts.get(Alteration.COPY_NUMBER);
		Change[] expz = alts.get(Alteration.EXPRESSION);

		if (muts == null || cncs == null || expz == null) return null;

		boolean[] x = new boolean[alts.getSize()];

		for (int i = 0; i < x.length; i++)
		{
			x[i] = (mutated ? muts[i].isAltered() : muts[i] == Change.NO_CHANGE) &&
				cncs[i] != Change.INHIBITING && expz[i] != Change.INHIBITING &&
				!cncs[i].isAbsent() && !expz[i].isAbsent();
		}

		return x;
	}

	public void calcDiffExpPvals()
	{
		targetChange = new HashMap<String, Map<String, Boolean>>();
		expPvals = new HashMap<String, Map<String, Double>>();
		expSigPvals = new HashMap<String, Map<String, Double>>();

		for (String mut : mutsig)
		{
			boolean[] normal = selectSubset(mut, false);
			boolean[] mutated = selectSubset(mut, true);

			if (normal == null || mutated == null) continue;

			Map<String, Double> pv = new HashMap<String, Double>();
			Map<String, Boolean> tc = new HashMap<String, Boolean>();

			for (String sym : graph.getSymbols())
			{
				double pval = calcDiffPval(sym, normal, mutated);

				if (Double.isNaN(pval)) continue;

				pv.put(sym, pval);
				tc.put(sym, Math.signum(calcMeanChange(sym, normal, mutated)) > 0);
			}

			expPvals.put(mut, pv);
			targetChange.put(mut, tc);

			List<String> result = FDR.select(pv, null, fdrThr);

			if (!result.isEmpty()) expSigPvals.put(mut, new HashMap<String, Double>());

			for (String gene : result)
			{
				expSigPvals.get(mut).put(gene, pv.get(gene));
			}
		}

		enrichSigPvals = new HashMap<String, Map<String, Double>>();

		for (String mut : expSigPvals.keySet())
		{
			enrichSigPvals.put(mut, new HashMap<String, Double>());
			Map<String, Double> pv = new HashMap<String, Double>();
			Map<String, Double> lim = new HashMap<String, Double>();

			for (String target : expSigPvals.get(mut).keySet())
			{
				Set<String> neighs = graph.getNeighbors(target);
				neighs = graph.getNeighbors(neighs);
				neighs.retainAll(expPvals.get(mut).keySet());

				int total = expPvals.get(mut).size();
				int featured = expSigPvals.get(mut).size();
				int selected = neighs.size();

				neighs.retainAll(expSigPvals.get(mut).keySet());

				int featuredSelected = neighs.size();

				pv.put(target, FishersExactTest.calcEnrichmentPval(
					total, featured, selected, featuredSelected));

				lim.put(target, FishersExactTest.calcEnrichmentPval(
					total, featured, selected, Math.min(selected, featured)));
			}

			List<String> select = FDR.select(pv, lim, fdrThr);

			for (String target : select)
			{
				enrichSigPvals.get(mut).put(target, pv.get(target));
			}
		}
	}

	public void calcComponents()
	{
		compsUp = new HashMap<String, List<Set<String>>>();
		compsDw = new HashMap<String, List<Set<String>>>();

		for (String mut : enrichSigPvals.keySet())
		{
			Set<String> up = new HashSet<String>();
			Set<String> dw = new HashSet<String>();

			for (String target : enrichSigPvals.get(mut).keySet())
			{
				if (targetChange.get(mut).get(target)) up.add(target);
				else dw.add(target);
			}

			ComponentSorter sorter = new ComponentSorter(up, graph);
			List<Set<String>> components = sorter.getComponents(componentSizeThr);
			if (!components.isEmpty()) compsUp.put(mut, components);

			sorter = new ComponentSorter(dw, graph);
			components = sorter.getComponents(componentSizeThr);
			if (!components.isEmpty()) compsDw.put(mut, components);
		}

		Set<String> resultGeneratingMuts = new HashSet<String>(compsUp.keySet());
		resultGeneratingMuts.addAll(compsDw.keySet());
		System.out.println("Result generated for " + resultGeneratingMuts.size() + " mutators");
	}

	public void writeComponents() throws IOException
	{
		writeComponents(compsUp, "-up");
		writeComponents(compsDw, "-dw");
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

					w2.write("node\t" + g1 + "\tcolor\t" + val2Color(expPvals.get(mut).get(g1),
						targetChange.get(mut).get(g1) ? 1 : -1) + "\n");
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
