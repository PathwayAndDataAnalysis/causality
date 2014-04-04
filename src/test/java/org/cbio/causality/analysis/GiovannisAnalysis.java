package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.binintanalysis.Dataset1;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.idmapping.HGNC;
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
public class GiovannisAnalysis
{
	private static final String dir = "binint/temp/";
	private Graph travSt;
	private Graph travExp;
	private Dataset1 dataset;
	private CBioPortalAccessor portalAcc;
	private ExpDataManager expMan;
	private Map<String, List<String>> downstream;
	private Map<String, List<String>> affectedDw;
	private int depth;
	double fdrThr;
	private Set<String> targetsToCheck;

	private static final boolean USE_COPY_NUMBER_UNCHANGED = false;
	private static final boolean ADD_CN_TO_MUT = false;
	private static final int MIN_GROUP_SIZE = 3;

	public static void main(String[] args) throws IOException
	{
		Graph graphSt;
		Graph graphExp;
		GiovannisAnalysis finder;
		Dataset1 dataset = Dataset1.THCA;
		int depth = 3;
		double fdrThr = 0.01;

		graphSt = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		graphExp = PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		graphExp.merge(MSigDBTFT.getGraph());
		finder = new GiovannisAnalysis(graphSt, graphExp, dataset, depth, fdrThr);

		checkZeros(finder.expMan);


		finder.find();
	}

	private static void checkZeros(ExpDataManager expMan)
	{
		int neverabove = 0;
		int sometimes = 0;
		int always = 0;
		Histogram h = new Histogram(0.1);
		for (String gene : HGNC.getAllSymbols())
		{
			if (!expMan.contains(gene)) continue;
			double r = expMan.getNonZeroRatio(gene);
			h.count(r);
			if (r == 0) neverabove++;
			else if (r == 1) always ++;
			else sometimes++;
		}
		System.out.println("neverabove = " + neverabove);
		System.out.println("always = " + always);
		System.out.println("sometimes = " + sometimes);
		h.print();
	}



	public List<String> find() throws IOException
	{
//		Set<String> g1 = new HashSet<String>(Arrays.asList("HRAS", "NRAS", "KRAS"));
//		Set<String> g2 = new HashSet<String>(Arrays.asList("BRAF"));

		Set<String> g1 = new HashSet<String>(Arrays.asList("HRAS"));
		Set<String> g2 = new HashSet<String>(Arrays.asList("NRAS", "KRAS"));

		Map<String, Double> pvals = calcDifferencePvals(g1, g2);

		System.out.println("pvals.size() = " + pvals.size());

		List<String> list = FDR.select(pvals, null, fdrThr);
		System.out.println("result size = " + list.size());

		Map<String, Boolean> dirs = getChangeDirections(list, g1, g2);

		int up = 0;
		int dw = 0;

		for (String gene : dirs.keySet())
		{
			Boolean d = dirs.get(gene);
			if (d) up++; else dw++;
		}

		System.out.println("up = " + up);
		System.out.println("dw = " + dw);

		return list;
	}

	private void decideTargets()
	{
		targetsToCheck = HGNC.getAllSymbols();
//		targetsToCheck = getDownstream("HRAS");
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

	public GiovannisAnalysis(Graph travSt, Graph travExp, Dataset1 dataset, int depth, double fdrThr)
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
		expMan.setTakeLog(true);
//		portalAcc.setCnVerifier(new CNVerifier(expMan, 0.05));
//		mutsig = BroadAccessor.getMutsigGenes(dataset.code(), 0.05);
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

	private double calcDiffPval(String symbol, boolean[] set1, boolean[] set2, boolean gaussian)
	{
		double[] exp = expMan.get(symbol);

//		if (exp == null) return Double.NaN;
//		if (Summary.variance(exp) < 2 || Summary.max(exp) < 3) return Double.NaN;

		boolean[] cnUnchanged = null;

		if (USE_COPY_NUMBER_UNCHANGED)
		{
			cnUnchanged = getCopyNumberUnchanged(symbol);
			if (cnUnchanged == null) return Double.NaN;
		}

		double[] vals1 = getSubset(exp, set1, cnUnchanged);
		double[] vals2 = getSubset(exp, set2, cnUnchanged);

		if (vals1.length < MIN_GROUP_SIZE || vals2.length < MIN_GROUP_SIZE)
		{
			return Double.NaN;
		}

		double pval;

		if (gaussian) pval = StudentsT.getPValOfMeanDifference(vals1, vals2);
		else pval = StudentsT.getPValOfMeanDifferenceBySimulation(vals1, vals2, 10000, 10);

		return pval;
	}

	private void printVals(double[] vals1, double[] vals2, String symbol)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter("/home/ozgun/Desktop/temp.txt"));

			writer.write("Ras-mutated\tBRAF-mutated\n");
			for (int i = 0; i < Math.max(vals1.length, vals2.length); i++)
			{
				if (i < vals1.length) writer.write(vals1[i] + "");
				writer.write("\t");
				if (i < vals2.length) writer.write(vals2[i] + "");
				writer.write("\n");
			}
			writer.close();

		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public Map<String, Boolean> getChangeDirections(List<String> genes, Set<String> g1, Set<String> g2)
	{
		boolean[][] b = getExclusiveMutations(g1, g2);
		Map<String, Boolean> dirs = new HashMap<String, Boolean>();
		for (String gene : genes)
		{
			dirs.put(gene, getChangeDirection(gene, b[0], b[1]));
		}
		return dirs;
	}

	private boolean getChangeDirection(String symbol, boolean[] set1, boolean[] set2)
	{
		double[] exp = expMan.get(symbol);

		if (exp == null) return false;

		boolean[] cnUnchanged = null;

		if (USE_COPY_NUMBER_UNCHANGED)
		{
			cnUnchanged = getCopyNumberUnchanged(symbol);
			if (cnUnchanged == null) return false;
		}

		double[] vals1 = getSubset(exp, set1, cnUnchanged);
		double[] vals2 = getSubset(exp, set2, cnUnchanged);

		return Summary.calcChangeOfMean(vals1, vals2) > 0;
	}

	private double[] getSubset(double[] vals, boolean[] inds, boolean[] cnUnchanged)
	{
		List<Double> list = new ArrayList<Double>(vals.length);

		for (int i = 0; i < vals.length; i++)
		{
			if (inds[i] && (cnUnchanged == null || cnUnchanged[i]) && !Double.isNaN(vals[i])) list.add(vals[i]);
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
		return getExpressionPvals(downstr, set1, set2);
	}

	private Map<String, Double> getExpressionPvals(Set<String> targets, boolean[] set1, boolean[] set2)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		Histogram2D h = new Histogram2D(0.2);


		System.out.println("targets = " + targets.size());
		Progress prg = new Progress(targets.size());
		for (String sym : targets)
		{
			prg.tick();

			if (expMan.getNonZeroRatio(sym) == 0) continue;
			double pval = calcDiffPval(sym, set1, set2, true);

			if (Double.isNaN(pval)) continue;
			if (pval == 0) pval = 1E-11;

			double pPerm = calcDiffPval(sym, set1, set2, false);
			if (pPerm == 0) pPerm = 1E-5;

			h.count(-Math.log(pval), -Math.log(pPerm));

			// pval = 0 is not real and it is not compatible with fisher's combined probability.
			// below is a better approximation.
//			if (pval == 0)
//			{
//				pval = 1E-11;
//			}

			map.put(sym, pval);
		}

		Histogram2DPlot p = new Histogram2DPlot(h);
		p.setLines(Arrays.asList(new double[]{1, 0}));
		p.setVisible(true);

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
			b[i] = !USE_COPY_NUMBER_UNCHANGED || !cnc[i].isAltered() && !cnc[i].isAbsent();
		}
		return b;
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



	public Map<String, Double> calcDifferingTargets(Set<String> g1, Set<String> g2)
	{
		boolean[][] b = getExclusiveMutations(g1, g2);

		return getExpressionPvals(targetsToCheck, b[0], b[1]);
	}

	private boolean[][] getExclusiveMutations(Set<String> g1, Set<String> g2)
	{
		boolean[][] b = new boolean[2][];

		b[0] = getMutated(g1);
		b[1] = getMutated(g2);

		System.out.println("b[0].length = " + b[0].length);
		System.out.println("M1 = " + ArrayUtil.countValue(b[0], true));
		System.out.println("M2 = " + ArrayUtil.countValue(b[1], true));

		for (int i = 0; i < b[0].length; i++)
		{
			if (b[0][i] && b[1][i])
			{
				b[0][i] = false;
				b[1][i] = false;
			}
		}

		System.out.println("After removing overlaps:");
		System.out.println("M1 = " + ArrayUtil.countValue(b[0], true));
		System.out.println("M2 = " + ArrayUtil.countValue(b[1], true));
		return b;
	}

	private boolean[] getMutated(Set<String> genes)
	{
		boolean[] b = new boolean[portalAcc.getAlterations(genes.iterator().next()).getSize()];

		for (String gene : genes)
		{
			AlterationPack alts = portalAcc.getAlterations(gene);
			Change[] ch = alts.get(Alteration.MUTATION);
			for (int i = 0; i < ch.length; i++)
			{
				if (ch[i].isAltered()) b[i] = true;
			}
		}
		return b;
	}

	public Map<String, Double> calcDifferencePvals(Set<String> g1, Set<String> g2)
	{
		decideTargets();
		return calcDifferingTargets(g1, g2);
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
