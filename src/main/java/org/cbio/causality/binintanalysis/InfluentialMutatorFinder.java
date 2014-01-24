package org.cbio.causality.binintanalysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.DownstreamTree;
import org.cbio.causality.analysis.GeneBranch;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.network.MSigDBTFT;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.util.*;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class InfluentialMutatorFinder
{
	private Graph travSt;
	private Graph travExp;
	private Dataset1 dataset;
	private CBioPortalAccessor portalAcc;
	private ExpDataManager expMan;
	private Map<String, List<String>> downstream;
	private Map<String, List<String>> affectedDw;
	private int depth;

	private static final int MIN_GROUP_SIZE = 1;

	public static void main(String[] args) throws IOException
	{
		InfluentialMutatorFinder finder = new InfluentialMutatorFinder(Dataset1.GBM, 2);

		double fdrThr = 0.05;

		Map<String, Double> pvals = finder.calcInfluencePvals();

		System.out.println("pvals.size() = " + pvals.size());

		List<String> list = FDR.select(pvals, null, fdrThr);
		System.out.println("result size = " + list.size());

		finder.generateInfluenceGraphs(list);

		for (String s : list)
		{
			System.out.println(s + "\t" + FormatUtil.roundToSignificantDigits(pvals.get(s), 2) +
				"\t" + finder.affectedDw.get(s));
		}
	}

	public InfluentialMutatorFinder(Dataset1 dataset, int depth)
		throws IOException
	{
		this.dataset = dataset;
		this.depth = depth;
		loadData();
		loadNetwork();
	}

	private void loadNetwork()
	{
//		travSt = SPIKE.getGraph();
		travSt = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
//		travSt.merge(SPIKE.getGraph());

//		travExp = MSigDBTFT.getGraph();
		travExp = PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		travExp.merge(MSigDBTFT.getGraph());
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

	private boolean getChangeDirection(String symbol, boolean[] set1, boolean[] set2)
	{
		double[] exp = expMan.get(symbol);

		if (exp == null) return false;

		boolean[] cnUnchanged = getCopyNumberUnchanged(symbol);
		if (cnUnchanged == null) return false;

		double[] vals1 = getSubset(exp, set1, cnUnchanged);
		double[] vals2 = getSubset(exp, set2, cnUnchanged);

		if (vals1.length < MIN_GROUP_SIZE || vals2.length < MIN_GROUP_SIZE) return false;

		return Summary.calcChangeOfMean(vals1, vals2) > 0;
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

	private Map<String, Double> getExpressionPvals(String mut, boolean[] set1, boolean[] set2)
	{
		return getExpressionPvals(getDownstream(mut), set1, set2);
	}

	private Map<String, Double> getExpressionPvals(Set<String> targets, boolean[] set1, boolean[] set2)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		for (String sym : targets)
		{
			double pval = calcDiffPval(sym, set1, set2);

			if (Double.isNaN(pval)) continue;

			map.put(sym, pval);
		}

		return map;
	}

	private Set<String> getDownstream(String mut)
	{
		Set<String> tfs = new HashSet<String>();
		tfs.add(mut);

		tfs.addAll(travSt.getDownstream(tfs, depth));

		return travExp.getDownstream(tfs);
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

	private boolean considerMutator(String symbol)
	{
		if (travExp.getDownstream(symbol).isEmpty()) return false;

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
		downstream = new HashMap<String, List<String>>();
		affectedDw = new HashMap<String, List<String>>();
		Map<String, Double> pvals = new HashMap<String, Double>();

		Set<String> symbols = travSt.getSymbols();
		symbols.addAll(travExp.getOneSideSymbols(true));

		for (String sym : symbols)
		{
			if (considerMutator(sym))
			{
				boolean[] normal = selectSubset(sym, false);
				boolean[] mutated = selectSubset(sym, true);

				final Map<String, Double> stat = getExpressionPvals(sym, normal, mutated);

				if (stat.size() < MIN_GROUP_SIZE) continue;

				List<String> list = new ArrayList<String>(stat.keySet());
				Collections.sort(list, new Comparator<String>()
				{
					@Override
					public int compare(String o1, String o2)
					{
						return stat.get(o1).compareTo(stat.get(o2));
					}
				});

				downstream.put(sym, list);

				list = FDR.select(stat, null, 0.05);
				affectedDw.put(sym, list);

				double pval = FishersCombinedProbability.pValue(toDoubleArray(stat));

				pvals.put(sym, pval);
			}
		}
		return pvals;
	}

	//--- Section: Graph generation----------------------------------------------------------------|

	private double calcPVal(String mut, Set<String> downs)
	{
		boolean[] normal = selectSubset(mut, false);
		boolean[] mutated = selectSubset(mut, true);

		Map<String, Double> pvals = getExpressionPvals(downs, normal, mutated);
		return FishersCombinedProbability.pValue(toDoubleArray(pvals));
	}

	private boolean calcChangeDirection(String mut, String target)
	{
		boolean[] normal = selectSubset(mut, false);
		boolean[] mutated = selectSubset(mut, true);

		return getChangeDirection(target, normal, mutated);
	}

	public void generateInfluenceGraphs(List<String> result)
	{
		DownstreamTree tree = new DownstreamTree(travSt, travExp);

		for (String mut : result)
		{
			GeneBranch g = tree.getTree(mut, downstream.get(mut), depth + 1);
			g.trimToMajorPaths(downstream.get(mut));
			if (g.selectLeaves(affectedDw.get(mut)))
			{
				writeTree(g, "temp");
			}
		}
	}

	private void writeTree(GeneBranch tree, String dir)
	{
		try
		{
			File d = new File(dir);
			if (!d.exists()) d.mkdirs();
			assert d.isDirectory();

			BufferedWriter writer1 = new BufferedWriter(new FileWriter(dir + "/" + tree.gene + ".sif"));
			BufferedWriter writer2 = new BufferedWriter(new FileWriter(dir + "/" + tree.gene + ".format"));

			writer2.write("graph\tgrouping\ton\n");
			writer2.write("node\t" + tree.gene + "\tcolor\t255 255 200\n");

			writeBranches(tree.gene, tree, writer1, writer2);

			writer1.close();
			writer2.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void writeBranches(String origTarget, GeneBranch branch, Writer writer1,
		Writer writer2) throws IOException
	{
		for (GeneBranch down : branch.branches)
		{
			writeSelf(origTarget, branch, down, writer1, writer2);
		}
	}

	private void writeSelf(String origTarget, GeneBranch currentParent, GeneBranch down,
		Writer writer1, Writer writer2) throws IOException
	{
		if (!down.isSelected()) return;

		String edgeTag;

		if (down.isLeaf() || travExp.getUpstream(down.gene).contains(currentParent))
		{
			edgeTag = SIFEnum.CONTROLS_EXPRESSION_OF.getTag();
		}
		else edgeTag = SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag();

		writer1.write(currentParent.gene + "\t" + edgeTag + "\t" + down.gene + "\n");

		writeWeights(origTarget, currentParent, down, edgeTag, writer2);

		writeBranches(origTarget, down, writer1, writer2);
	}

	private void writeWeights(String orig, GeneBranch from, GeneBranch to, String edgeType,
		Writer writer) throws IOException
	{
		Set<String> dwstr = to.getAllGenes();
		dwstr.retainAll(downstream.get(orig));
		assert !dwstr.isEmpty();
		double cumPval = calcPVal(orig, dwstr);
		boolean upreg = calcChangeDirection(orig, to.gene);


		String key = from.gene + " " + edgeType + " " + to.gene;
		writer.write("edge\t" + key + "\tcolor\t" + val2Color(cumPval, 0) + "\n");
		writer.write("edge\t" + key + "\twidth\t2\n");

		if (affectedDw.get(orig).contains(to.gene))
		{
			double pval = calcPVal(orig, Collections.singleton(to.gene));
			writer.write("node\t" + to.gene + "\tcolor\t" + val2Color(pval, upreg ? 1 : -1) + "\n");
		}
		else
		{
			writer.write("node\t" + to.gene + "\tcolor\t255 255 255\n");
		}
	}

	private String val2Color(double pval, int type)
	{
		double score = Math.min(10, -Math.log(pval));

		int v = 255 - (int) Math.round((255D / 10D) * score);

		switch (type)
		{
			case 0: return  v + " " + v + " " + v;
			case 1: return "255 " + v + " " + v;
			case -1: return v + " 255 " + v;
			default: throw new IllegalArgumentException();
		}
	}

}
