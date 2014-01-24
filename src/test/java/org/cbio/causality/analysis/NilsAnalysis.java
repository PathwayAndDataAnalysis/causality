package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersExactTest;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class NilsAnalysis
{
	@Test
	@Ignore
	public void printDegreeDistribution() throws IOException
	{
		Graph trav = new Graph();
		trav.load("/home/ozgun/Projects/biopax-pattern/Related-through-interaction.txt",
			new HashSet<String>(Arrays.asList(SIFEnum.NEIGHBOR_OF.getTag())),
			Collections.<String>emptySet());

		Set<String> pancan = readGenes("pancan12_drivers_not_in_variants");
		Set<String> variant = readGenes("variant_list_filtered.genes");
		System.out.println("variant.size() = " + variant.size());

		Set<String> symbols = trav.getSymbols();
		variant.retainAll(symbols);
		System.out.println("variant.size() = " + variant.size());

		final Map<String, Integer> hits = new HashMap<String, Integer>();
		final Map<String, Double> pvals = new HashMap<String, Double>();
		final Map<String, Double> limits = new HashMap<String, Double>();

		for (String panGene : pancan)
		{
			Set<String> neighbors = trav.getNeighbors(panGene);
			int score = intersectCnt(neighbors, variant);

			hits.put(panGene, score);

			Set<String> n = new HashSet<String>(symbols);
			n.remove(panGene);

			Set<String> s00 = new HashSet<String>(n);
			s00.removeAll(neighbors);
			s00.removeAll(variant);

			Set<String> s01 = new HashSet<String>(n);
			s01.removeAll(neighbors);
			s01.retainAll(variant);

			Set<String> s10 = new HashSet<String>(n);
			s10.retainAll(neighbors);
			s10.removeAll(variant);

			Set<String> s11 = new HashSet<String>(n);
			s11.retainAll(neighbors);
			s11.retainAll(variant);

			int n00 = s00.size();
			int n10 = s10.size();
			int n01 = s01.size();
			int n11 = s11.size();

			double pv = FishersExactTest.calcPositiveDepPval(n00, n01, n10, n11);

			pvals.put(panGene, pv);

			n11 += n10;
			n00 += n10;
			n01 -= n10;
			assert n01 >= 0;
			n10 = 0;

			pv = FishersExactTest.calcPositiveDepPval(n00, n01, n10, n11);
			limits.put(panGene, pv);
		}

		Map<String, Double> qvals = FDR.getQVals(pvals, limits);

		List<String> sorted = new ArrayList<String>(pancan);

		Collections.sort(sorted, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return pvals.get(o1).compareTo(pvals.get(o2));
			}
		});

		DecimalFormat f = new DecimalFormat("0.000");
		BufferedWriter writer = new BufferedWriter(new FileWriter("NilsRankedList.txt"));
		writer.write("Gene\tHit\tP-value\tQ-value");

		for (String gene : sorted)
		{
			writer.write("\n" + gene + "\t" + hits.get(gene) +
				"\t" + f.format(pvals.get(gene)) +
				"\t" + f.format(qvals.get(gene)));
		}

		writer.close();
	}

	private int intersectCnt(Set<String> small, Set<String> large)
	{
		int cnt = 0 ;

		for (String s : small)
		{
			if (large.contains(s)) cnt++;
		}
		return cnt;
	}

	private Set<String> readGenes(String file)
	{
		Scanner sc = new Scanner(
			getClass().getResourceAsStream(file));

		Set<String> genes = new HashSet<String>();

		while(sc.hasNextLine())
		{
			String gene = sc.nextLine();

			if (!gene.isEmpty()) genes.add(gene);
		}
		return genes;
	}
}
