package org.cbio.causality.binintanalysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.*;
import org.cbio.causality.data.go.GO;
import org.cbio.causality.data.portal.BroadAccessor;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.network.*;
import org.cbio.causality.util.*;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class InfluentialMutatorFinder
{
	private static final Set<String> focusExp = new HashSet<String>(Arrays.asList(
		"ERBB4", "GRIN2A", "GRIN2B", "RELN", "DAB1"));

	private static final Set<String> focusMut = null;
//	private static final Set<String> focusMut = new HashSet<String>(Arrays.asList(
//		"TP53", "BRAF", "NRAS", "UGT2B15", "TTN", "CDKN2A"));

	private static final String dir = "binint/InfluentialMutatorFinder/";
	private Graph travSt;
	private Graph travExp;
	private Dataset1 dataset;
	private CBioPortalAccessor portalAcc;
	private ExpDataManager expMan;
	private Map<String, List<String>> downstream;
	private Map<String, List<String>> affectedDw;
	private int depth;
	Set<String> mutsig;
	double fdrThr;

	private static final boolean LIMIT_TO_MUTSIG = true;
	private static final boolean ADD_CN_TO_MUT = false;
	private static final int MIN_GROUP_SIZE = 2;
	private static final double MIN_MUT_RATIO = 0.0;

	public static void main(String[] args) throws IOException
	{
		Graph graphSt;
		Graph graphExp;
		InfluentialMutatorFinder finder;
		Dataset1 dataset = Dataset1.SKCM;
		int depth = 3;
		double fdrThr = 0.05;

		graphSt = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		graphExp = PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		graphExp.merge(MSigDBTFT.getGraph());
		finder = new InfluentialMutatorFinder(graphSt, graphExp, dataset, depth, fdrThr);
		List<String> resPC = finder.find();

		if (true) return;

		graphSt = SPIKE.getGraphPostTl();
		graphExp = SPIKE.getGraphTR();
		graphExp.merge(MSigDBTFT.getGraph());
		finder = new InfluentialMutatorFinder(graphSt, graphExp, dataset, depth, fdrThr);
		List<String> resSpike = finder.find();

		graphSt = SignaLink.getGraphPostTl();
		graphExp = SignaLink.getGraphTR();
		graphExp.merge(MSigDBTFT.getGraph());
		finder = new InfluentialMutatorFinder(graphSt, graphExp, dataset, depth, fdrThr);
		List<String> resSignalink = finder.find();

		graphSt = PathwayCommons.getGraph(SIFEnum.INTERACTS_WITH).changeTo(true);
		graphSt.merge(HPRD.getGraph(true));
		graphSt.merge(IntAct.getGraph(true));
		graphExp = MSigDBTFT.getGraph();
		finder = new InfluentialMutatorFinder(graphSt, graphExp, dataset, depth, fdrThr);
		List<String> resPPI = finder.find();

//		graphSt.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
//		graphSt.merge(SPIKE.getGraphPostTl());
//		graphExp.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
//		graphExp.merge(SPIKE.getGraphTR());
//		graphExp.merge(MSigDBTFT.getGraph());
//		finder = new InfluentialMutatorFinder(graphSt, graphExp, dataset, depth, fdrThr);
//		List<String> merged = finder.find();

		CollectionUtil.printVennCounts(resPC, resSpike, resSignalink, resPPI);
	}
	public List<String> find() throws IOException
	{
		Map<String, Double> pvals = calcInfluencePvals();

		System.out.println("pvals.size() = " + pvals.size());

		List<String> list = FDR.select(pvals, null, fdrThr);
		System.out.println("result size = " + list.size());

		generateInfluenceGraphs(list);

		List<String> mutsig = compareMutSig(list, true);
		List<String> nonmut = compareMutSig(list, false);

		System.out.println("\nMutsig genes (" + mutsig.size() + " of " + this.mutsig.size() + ")");
		printResults(pvals, mutsig);
		System.out.println("\nNon mutsig genes (" + nonmut.size() + ")");
		printResults(pvals, nonmut);

		GO.printEnrichment(new HashSet<String>(list), pvals.keySet(), 0.05);
		return list;
	}

	private void printResults(Map<String, Double> pvals,
		List<String> list)
	{
		for (String s : list)
		{
			System.out.println(s + "\t" + FormatUtil.roundToSignificantDigits(pvals.get(s), 2) +
				"\t" + affectedDw.get(s));
		}
	}

	public InfluentialMutatorFinder(Graph travSt, Graph travExp, Dataset1 dataset, int depth, double fdrThr)
		throws IOException
	{
		this.travSt = travSt;
		this.travExp = travExp;
		this.dataset = dataset;
		this.depth = depth;
		this.fdrThr = fdrThr;
		loadData();
	}

	private void loadData() throws IOException
	{
		portalAcc = new CBioPortalAccessor(dataset.mutCnCallExpZ);
		expMan = new ExpDataManager(portalAcc.getGeneticProfileById(dataset.exp.profileID[0]),
			portalAcc.getCaseListById(dataset.exp.caseListID));
		portalAcc.setCnVerifier(new CNVerifier(expMan, 0.05));

		mutsig = BroadAccessor.getMutsigGenes(dataset.code(), 0.05);
	}

	private List<String> compareMutSig(List<String> list, boolean intersect)
	{
		List<String> select = new ArrayList<String>(list);
		if (intersect) select.retainAll(mutsig); else select.removeAll(mutsig);
		return select;
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
			if (ADD_CN_TO_MUT)
			{
				x[i] = (mutated ? muts[i].isAltered() || cncs[i].isAltered() :
					muts[i] == Change.NO_CHANGE) && cncs[i] == Change.NO_CHANGE && expz[i] !=
					Change.INHIBITING && !expz[i].isAbsent();
			}
			else
			{
				x[i] = (mutated ? muts[i].isAltered() : muts[i] == Change.NO_CHANGE) &&
					cncs[i] != Change.INHIBITING && expz[i] != Change.INHIBITING &&
					!cncs[i].isAbsent() && !expz[i].isAbsent();
			}
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
		Set<String> downstr = getDownstream(mut);

		if (focusExp != null) downstr.retainAll(focusExp);

		return getExpressionPvals(downstr, set1, set2);
	}

	private Map<String, Double> getExpressionPvals(Set<String> targets, boolean[] set1, boolean[] set2)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		for (String sym : targets)
		{
			double pval = calcDiffPval(sym, set1, set2);

			if (Double.isNaN(pval)) continue;

			// pval = 0 is not real and it is not compatible with fisher's combined probability.
			// below is a better approximation.
			if (pval == 0)
			{
				pval = 1E-11;
			}

			map.put(sym, pval);
		}

		return map;
	}

	private Set<String> getDownstream(String mut)
	{
		Set<String> tfs = new HashSet<String>();
		tfs.add(mut);

		if (depth > 1) tfs.addAll(travSt.getDownstream(tfs, depth - 1));

		return travExp.getDownstream(tfs);
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

	private boolean considerMutator(String symbol)
	{
		if (getDownstream(symbol).isEmpty()) return false;

		AlterationPack alts = portalAcc.getAlterations(symbol);

		if (alts == null) return false;

		if (alts.countAltered(Alteration.MUTATION) < MIN_GROUP_SIZE) return false;
		if (alts.getAlteredRatio(Alteration.MUTATION) < MIN_MUT_RATIO) return false;
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

		if (focusMut != null) symbols.retainAll(focusMut);
		if (LIMIT_TO_MUTSIG) symbols.retainAll(mutsig);

		for (String sym : symbols)
		{
			if (considerMutator(sym))
			{
				boolean[] normal = selectSubset(sym, false);
				boolean[] mutated = selectSubset(sym, true);

				final Map<String, Double> stat = getExpressionPvals(sym, normal, mutated);

				if (stat.isEmpty()) continue;

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
		DownstreamTree tree = new DownstreamTree(travSt, travExp, new BranchDataProvider()
		{
			@Override
			public Color getColor(String gene, String root)
			{
				double pval = calcPVal(root, Collections.singleton(gene));

				String colS = val2Color(pval, calcChangeDirection(root, gene) ? 1 : -1);
				String[] c = colS.split(" ");
				return new Color(Integer.parseInt(c[0]), Integer.parseInt(c[1]),
					Integer.parseInt(c[2]));
			}

			@Override
			public double getThickness(GeneBranch branch, String root)
			{
				Set<String> dwstr = branch.getAllGenes();
				double cumPval = calcPVal(root, dwstr);

				return Math.min(-Math.log(cumPval), 5);
			}
		});

		String dir  = this.dir + dataset.name() + "/";
		File d = new File(dir);
		if (!d.exists()) d.mkdirs();
		assert d.isDirectory();

		for (String mut : result)
		{
			GeneBranch g = tree.getTree(mut, downstream.get(mut), depth);
			g.trimToMajorPaths(downstream.get(mut));
			if (g.selectLeaves(affectedDw.get(mut)))
			{
				if (!affectedDw.get(mut).isEmpty() && affectedDw.get(mut).size() < 6)
					RadialInfluenceTree.write(g.copy(true), false, dir + mut + ".svg");

				writeTree(g, dir);
			}
		}
	}

	private void writeTree(GeneBranch tree, String dir)
	{
		try
		{
			BufferedWriter writer1 = new BufferedWriter(new FileWriter(dir + tree.gene + ".sif"));
			BufferedWriter writer2 = new BufferedWriter(new FileWriter(dir + tree.gene + ".format"));

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
