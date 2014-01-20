package org.cbio.causality.binintanalysis;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.*;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class InfluentialTFFinder
{
	private Graph trav;
	private Dataset1 dataset;
	CBioPortalAccessor portalAcc;
	ExpDataManager expMan;
	Map<String, List<String>> affectedDw;

	private static final int MIN_GROUP_SIZE = 1;

	public static void main(String[] args) throws IOException
	{
		InfluentialTFFinder finder = new InfluentialTFFinder(Dataset1.THCA,
			"/home/ozgun/Projects/chibe/portal-cache/PC.sif");

		double fdrThr = 0.05;

		Map<String, Double> result = finder.calcInfluencePvals();

		System.out.println("pvals.size() = " + result.size());

		List<String> list = FDR.select(result, null, fdrThr);
		System.out.println("result size = " + list.size());

		for (String s : list)
		{
			System.out.println(s + "\t" + result.get(s) + "\t" + finder.affectedDw.get(s));
		}
	}

	public InfluentialTFFinder(Dataset1 dataset, String networkFile)
		throws IOException
	{
		this.dataset = dataset;
		loadData();
		loadNetwork(networkFile);
	}

	private void loadNetwork(String filename)
	{
		trav = new Graph();
		trav.load(filename, Collections.<String>emptySet(),
			Collections.singleton(SIFType.CONTROLS_EXPRESSION_OF.getTag()));
	}

	private void loadData() throws IOException
	{
		portalAcc = new CBioPortalAccessor(dataset.mutCnCallExpZ);
		expMan = new ExpDataManager(portalAcc.getGeneticProfileById(dataset.exp.profileID[0]),
			portalAcc.getCaseListById(dataset.exp.caseListID));
	}

	/**
	 * Gets the samples that the gene is either mutated or intact. In both cases it is expressed,
	 * and not copy number lost.
	 */
	private boolean[] selectSubset(String symbol, boolean mutated)
	{
		AlterationPack alts = portalAcc.getAlterations(symbol);

		Change[] muts = alts.get(Alteration.MUTATION);
		Change[] cncs = alts.get(Alteration.COPY_NUMBER);
		Change[] expz = alts.get(Alteration.EXPRESSION);

		boolean[] x = new boolean[alts.getSize()];

		for (int i = 0; i < x.length; i++)
		{
			x[i] = (mutated ? muts[i].isAltered() : muts[i] == Change.NO_CHANGE) &&
				cncs[i] != Change.INHIBITING && expz[i] != Change.INHIBITING &&
				!cncs[i].isAbsent() && !expz[i].isAbsent();
		}

		return x;
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

	private Map<String, Double> getExpressionPvals(String tf, boolean[] set1, boolean[] set2)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		for (String sym : trav.getDownstream(tf))
		{
			double pval = calcDiffPval(sym, set1, set2);

			if (Double.isNaN(pval)) continue;

			map.put(sym, pval);
		}

		return map;
	}

	private boolean[] getCopyNumberUnchanged(String symbol)
	{
		AlterationPack alts = portalAcc.getAlterations(symbol);
		Change[] cnc = alts.get(Alteration.COPY_NUMBER);

		if (cnc == null) return null;

		boolean[] b = new boolean[cnc.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !cnc[i].isAltered() && !cnc[i].isAbsent();
		}
		return b;
	}

	private boolean considerTF(String symbol)
	{
		if (trav.getDownstream(symbol).isEmpty()) return false;

		AlterationPack alts = portalAcc.getAlterations(symbol);

		if (alts == null) return false;

		if (alts.countAltered(Alteration.MUTATION) < MIN_GROUP_SIZE) return false;
		if (alts.isAbsent(Alteration.COPY_NUMBER)) return false;
		if (alts.isAbsent(Alteration.EXPRESSION)) return false;

		return true;
	}

	private double[] toDoubleArray(Map<String, Double> map)
	{
		double[] d = new double[map.size()];

		int i = 0;
		for (String key : map.keySet())
		{
			d[i++] = map.get(key);
		}
		return d;
	}

	public Map<String, Double> calcInfluencePvals()
	{
		affectedDw = new HashMap<String, List<String>>();
		Map<String, Double> pvals = new HashMap<String, Double>();

		for (String sym : trav.getSymbols())
		{
			if (considerTF(sym))
			{
				boolean[] normal = selectSubset(sym, false);
				boolean[] mutated = selectSubset(sym, true);

				final Map<String, Double> stat = getExpressionPvals(sym, normal, mutated);

				if (stat.size() < MIN_GROUP_SIZE) continue;

				ArrayList<String> list = new ArrayList<String>(stat.keySet());
				Collections.sort(list, new Comparator<String>()
				{
					@Override
					public int compare(String o1, String o2)
					{
						return stat.get(o1).compareTo(stat.get(o2));
					}
				});
				affectedDw.put(sym, list);

				double pval = FishersCombinedProbability.pValue(toDoubleArray(stat));

				pvals.put(sym, pval);
			}
		}
		return pvals;
	}
}
