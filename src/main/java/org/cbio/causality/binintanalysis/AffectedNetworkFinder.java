package org.cbio.causality.binintanalysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.BroadAccessor;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;
import org.cbio.causality.util.FormatUtil;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class AffectedNetworkFinder
{
	private Graph trav;
	private int totalSymbolSize;
	private String study;
	private double fdrThr;
	Set<String> mutsig;

	public static void main(String[] args) throws IOException
	{
		AffectedNetworkFinder anf = new AffectedNetworkFinder("COADREAD", 0.05, 0.25);

		Map<String, Double> pvals = anf.calcPvals();

		anf.writeNetwork(pvals, "/home/ozgun/Desktop/temp.sif");
	}

	public AffectedNetworkFinder(String study, double mutsigThr, double fdrThr)
	{
		this.study = study;
		loadNetwork();
		mutsig = BroadAccessor.getMutsigGenes(study, mutsigThr);
		mutsig.retainAll(trav.getSymbols());
		System.out.println("mutsig in network = " + mutsig.size());
		this.fdrThr = fdrThr;
	}

	public Map<String, Double> calcPvals()
	{
		Map<String, Double> scores = new HashMap<String, Double>();

		for (String sym : trav.getSymbols())
		{
			double pv = calcUpstreamEnrichmentPval(sym, 1);
			scores.put(sym, pv);
		}
		return scores;
	}

	private double calcUpstreamEnrichmentPval(String symbol, int depth)
	{
		double minPval = 1;
		for (int i = 0; i < depth; i++)
		{
			Set<String> up = collectUpstream(symbol, i + 1);

			int cnt = 0;

			for (String s : up)
			{
				if (mutsig.contains(s)) cnt++;
			}

			double pval = FishersExactTest.calcEnrichmentPval(
				totalSymbolSize, mutsig.size(), up.size(), cnt);

			if (pval < minPval) minPval = pval;
		}
		return minPval;
	}

	private Set<String> collectUpstream(String sym, int depth)
	{
		Set<String> up = new HashSet<String>(trav.getUpstream(sym));
		Set<String> newG = new HashSet<String>(up);

		for (int i = 1; i < depth; i++)
		{
			Set<String> newUp = trav.getUpstream(newG);
			newUp.removeAll(newG);

			newG.clear();
			newG.addAll(newUp);
			newG.removeAll(up);
			
			up.addAll(newUp);			
		}
		return up;
	}

	private void loadNetwork()
	{
		trav = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		trav.merge(SPIKE.getGraph());
		totalSymbolSize = trav.getSymbols().size();
	}

	private void writeNetwork(Map<String, Double> pvals, String outFile) throws IOException
	{
		List<String> selected = new ArrayList<String>(FDR.select(pvals, null, fdrThr));

		System.out.println("selected.size() = " + selected.size());
		for (String gene : selected)
		{
			System.out.println(gene + "\t" +
				FormatUtil.roundToSignificantDigits(pvals.get(gene), 2));
		}

		selected.addAll(mutsig);
		writeNetwork(selected, outFile, pvals);
	}

	private void writeNetwork(Collection<String> selected, String filename, Map<String, Double> pvals) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (String s : selected)
		{
			for (String dw : trav.getDownstream(s))
			{
				if (selected.contains(dw))
				{
					writer.write(s + "\t" + SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag() + "\t" + dw +
						"\n");
				}
			}
		}

		writer.close();

//		double maxScore = maxScore(pvals);

		writer = new BufferedWriter(new FileWriter(filename.substring(0, filename.lastIndexOf(".")) + ".format"));

		for (String s : selected)
		{
			if (!pvals.containsKey(s)) continue;

//			int v = getColorOfScore(scores.get(s), maxScore);
//			writer.write(s + "\tcolor\t255 " + v + " " + v + "\n");
			writer.write("node\t" + s + "\tcolor\t255 255 255\n");
			if (mutsig.contains(s))
				writer.write("node\t" + s + "\thighlight\ton\n");
		}

		writer.close();

	}

	private int getColorOfScore(double val, double max)
	{
		return 255 - (int) Math.round((val / max) * 255);
	}

	private double maxScore(Map<String, Double> scores)
	{
		double max = -1;
		for (Double score : scores.values())
		{
			if (score > max) max = score;
		}
		return max;
	}
}
