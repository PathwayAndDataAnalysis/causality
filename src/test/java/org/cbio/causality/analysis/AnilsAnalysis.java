package org.cbio.causality.analysis;

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
public class AnilsAnalysis
{
	@Test
	@Ignore
	public void printRankedNeighborhood() throws IOException
	{
		Traverse trav = new Traverse();
		trav.load("/home/ozgun/Projects/biopax-pattern/Related-through-interaction.txt",
			new HashSet<String>(Arrays.asList(SIFType.RELATED_THROUGH_INTERACTION.getTag())),
			Collections.<String>emptySet());

		Set<String> seeds = readGenes("resistance-genes.txt");

		Set<String> neighs = trav.getNeighbors(seeds);
		neighs.removeAll(seeds);

		final Map<String, Integer> scores = new HashMap<String, Integer>();

		for (String neigh : neighs)
		{
			int cnt = 0;

			for (String seed : seeds)
			{
				if (trav.getNeighbors(seed).contains(neigh)) cnt++;
			}

			scores.put(neigh, cnt);
		}

		List<String> sorted = new ArrayList<String>(neighs);

		Collections.sort(sorted, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return scores.get(o2).compareTo(scores.get(o1));
			}
		});

		BufferedWriter writer = new BufferedWriter(new FileWriter("ranked-neighborhood.txt"));
		writer.write("Score\tGene");

		for (String neigh : sorted)
		{
			writer.write("\n" + scores.get(neigh) + "\t" + neigh);
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
