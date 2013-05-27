package org.cbio.causality.alteredreaction;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CBioPortalManager;
import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.hprd.HPRD;
import org.cbio.causality.hprd.InteractionProvider;
import org.cbio.causality.hprd.ShuffledHPRD;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Progress;
import org.cbio.causality.util.ScoreUtil;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class CoocMutAnalyzer
{
	CBioPortalManager cman;
	CBioPortalAccessor cacc;
	List<Reaction> reactions;

	Map<String, Double> reacScore;
	Map<Interaction, Interaction> interScore;
	Map<String, AlterationPack> cache;

	Set<String> genes;

	public CoocMutAnalyzer(PortalDataset dataset)
	{
		try{cacc = new CBioPortalAccessor(dataset);
		} catch (IOException e){e.printStackTrace();}

		cman = new CBioPortalManager();

		AlteredReactionFinder arf = new AlteredReactionFinder();
		reactions = arf.readReactions();

		reacScore = new HashMap<String, Double>();
		cache = new HashMap<String, AlterationPack>();
	}

	public void resetScores()
	{
		if (reacScore != null) reacScore.clear();
		if (interScore != null) interScore.clear();
	}

	public AlterationPack getAlterationPack(String symbol)
	{
		if (cache.containsKey(symbol)) return cache.get(symbol);

		String[] data = cman.getDataForGene(symbol,
			cacc.getCurrentGeneticProfiles().iterator().next(), cacc.getCurrentCaseList());

		if (data == null)
		{
			cache.put(symbol, null);
			return null;
		}

		Change[] ch = new Change[data.length];
		for (int i = 0; i < data.length; i++)
		{
			if (!data[i].isEmpty() && !data[i].equals("NaN"))
			{
				ch[i] = Change.UNKNOWN_CHANGE;
			}
			else
			{
				ch[i] = Change.NO_CHANGE;
			}
		}

		AlterationPack pack = new AlterationPack(symbol);
		pack.put(Alteration.MUTATION, ch);

		cache.put(symbol, pack);
		return pack;
	}

	public void scoreReactionCoMutations()
	{
		for (Reaction r : reactions)
		{
			List<AlterationPack> genes = new ArrayList<AlterationPack>();

			for (String gene : r.genes)
			{
				AlterationPack pack = getAlterationPack(gene);
				if (pack != null && pack.isAltered()) genes.add(pack);
			}

			if (genes.size() < 2)
			{
				reacScore.put(r.ID, 0D);
				continue;
			}

			int score = 0;
			for (int i = 0; i < genes.get(0).getSize(); i++)
			{
				int cnt = 0;
				for (AlterationPack gene : genes)
				{
					if (gene.getChange(Alteration.MUTATION, i).isAltered()) cnt++;
				}
				if (cnt > 1) score++;
			}

			reacScore.put(r.ID, score / (double) r.genes.size());
		}
	}

	public void sortReacToScore()
	{
		Collections.sort(reactions, new Comparator<Reaction>()
		{
			@Override
			public int compare(Reaction o1, Reaction o2)
			{
				return reacScore.get(o2.ID).compareTo(reacScore.get(o1.ID));
			}
		});
	}

	public void writeReactionsWithScore(String filename) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (Reaction r : reactions)
		{
			writer.write("\n" + reacScore.get(r.ID) + "\t" + r);
		}

		writer.close();
	}

	public void shuffleGeneIdentitiesInReac()
	{
		Set<String> names = new HashSet<String>();
		for (Reaction r : reactions)
		{
			names.addAll(r.genes);
		}

		List<String> list1 = new ArrayList<String>(names);
		List<String> list2 = new ArrayList<String>(names);

		Collections.sort(list1);
		Collections.shuffle(list2);

		Map<String, String> map = new HashMap<String, String>();

		for (int i = 0; i < names.size(); i++)
		{
			map.put(list1.get(i), list2.get(i));
		}

		for (Reaction r : reactions)
		{
			ArrayList<String> genes = new ArrayList<String>();
			for (String g : r.genes)
			{
				genes.add(map.get(g));
			}
			r.genes = genes;
		}
	}

	public Histogram getReacScoreDistribution()
	{
		Histogram h = new Histogram(1);
		for (Double d : reacScore.values())
		{
			h.count(d);
		}
		return h;
	}

	public void prepareInteractingGenes()
	{
		interScore = new HashMap<Interaction, Interaction>();

		Set<String> all = HPRD.getAllSymbols();
		genes = new HashSet<String>();
		Set<String> seq = new HashSet<String>();

		for (String gene : all)
		{
			AlterationPack pack = getAlterationPack(gene);
			if (pack != null)
			{
				seq.add(gene);
				if (pack.isAltered()) genes.add(gene);
			}
		}
		HPRD.cropTo(genes);

		System.out.println("seq.size() = " + seq.size());
		System.out.println("altered.size() = " + genes.size());
	}

	public void scoreInteractionCoMutations(InteractionProvider ip)
	{
		int size = getAlterationPack(genes.iterator().next()).getSize();

		for (int i = 0; i < size; i++)
		{
			Set<String> set = new HashSet<String>();

			for (String gene : genes)
			{
				if (getAlterationPack(gene).getChange(Alteration.MUTATION, i).isAltered())
				{
					set.add(gene);
				}
			}

			if (set.size() > 1)
			{
				for (String s1 : set)
				{
					Set<String> inter = ip.getInteractions(s1);
					inter.retainAll(set);

					for (String s2 : inter)
					{
						if (s1.compareTo(s2) >= 0) continue;

						Interaction in = new Interaction(s1, s2);

						if (!interScore.containsKey(in)) interScore.put(in, in);
						interScore.get(in).score++;
					}
				}
			}
		}
	}

	public Histogram getInterScoreHisto()
	{
		Histogram h = new Histogram(1);
		for (Interaction in : interScore.keySet())
		{
			h.count(in.score);
		}
		return h;
	}

	public ScoreUtil getInterScores()
	{
		ScoreUtil su = new ScoreUtil();
		for (Interaction in : interScore.keySet())
		{
			su.addSCore(in.score);
		}
		return su;
	}

	public void writeCoocInterAsSIF(String filename, double thr) throws IOException
	{
		List<Interaction> inter = new ArrayList<Interaction>(interScore.keySet());
		Collections.sort(inter);
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		int i = 0;
		for (Interaction in : inter)
		{
			if (in.score >= thr)
			{
				i++;
				writer.write(in.sym1 + "\tBINDS_TO\t" + in.sym2 + "\n");
			}
		}

		writer.close();
		System.out.println("result interactions = " + i);
	}

	public void calcNeighborhoodTendency(InteractionProvider ip)
	{
		Map<String, Integer> coMutCnt = new HashMap<String, Integer>();
		Map<String, Integer> neighCnt = new HashMap<String, Integer>();

		int size = getAlterationPack(genes.iterator().next()).getSize();

		for (int i = 0; i < size; i++)
		{
			Set<String> set = new HashSet<String>();

			for (String gene : genes)
			{
				if (getAlterationPack(gene).getChange(Alteration.MUTATION, i).isAltered())
				{
					set.add(gene);
				}
			}

			if (set.size() > 1)
			{
				for (String s1 : set)
				{
					Set<String> inter = ip.getInteractions(s1);

					if (!coMutCnt.containsKey(s1)) coMutCnt.put(s1, 0);
					coMutCnt.put(s1,  coMutCnt.get(s1) + set.size() - 1);

					inter.retainAll(set);

					if (!neighCnt.containsKey(s1)) neighCnt.put(s1, 0);
					neighCnt.put(s1,  neighCnt.get(s1) + inter.size());
				}
			}
		}

		int cnt1 = 0;
		int cnt2 = 0;

		for (String gene : genes)
		{
			int comut = coMutCnt.get(gene);
			int neighComut = neighCnt.get(gene);
			int neigh = ip.getInteractions(gene).size();
			int total = genes.size() - 1;

			double nrat = neigh / (double) total;
			double mrat = neighComut / (double) comut;

			if (mrat > nrat) cnt1++;
			else cnt2++;
		}

		System.out.println("cnt1 = " + cnt1);
		System.out.println("cnt2 = " + cnt2);
	}

	public double[] getUnexpectedDistribution(InteractionProvider ip)
	{
		double[] s = new double[2];
		for (String gene1 : genes)
		{
			AlterationPack pack1 = getAlterationPack(gene1);

			for (String gene2 : ip.getInteractions(gene1))
			{
				if (!genes.contains(gene2) || gene1.compareTo(gene2) >= 0) continue;

				AlterationPack pack2 = getAlterationPack(gene2);

				int[] cnt = getCounts(pack1.get(Alteration.MUTATION), pack2.get(Alteration.MUTATION));

				double e = (cnt[1] * cnt[2]) / (double) pack1.getSize();

				double dif = cnt[0] - e;

				if (dif > 0) s[0] += dif;
				else if (dif < 0) s[1] -= dif;
				else System.out.println("Equal!!");
			}
		}
		return s;
	}

	public double[] getAllPairsDistribution()
	{
		double[] s = new double[2];
		for (String gene1 : genes)
		{
			AlterationPack pack1 = getAlterationPack(gene1);

			for (String gene2 : genes)
			{
				if (gene1.compareTo(gene2) >= 0) continue;

				AlterationPack pack2 = getAlterationPack(gene2);

				int[] cnt = getCounts(pack1.get(Alteration.MUTATION), pack2.get(Alteration.MUTATION));

				double e = (cnt[1] * cnt[2]) / (double) pack1.getSize();

				double dif = cnt[0] - e;

				if (dif > 0) s[0] += dif;
				else if (dif < 0) s[1] -= dif;
				else System.out.println("Equal!!");
			}
		}
		return s;
	}

	private int[] getCounts(Change[] ch1, Change[] ch2)
	{
		int[] c = new int[3];

		for (int i = 0; i < ch1.length; i++)
		{
			if (ch1[i].isAltered())
			{
				c[1]++;
				if (ch2[i].isAltered())
				{
					c[2]++;
					c[0]++;
				}
			}
			else if (ch2[i].isAltered()) c[2]++;
		}
		return c;
	}

	public static void processReac(String[] args) throws IOException
	{
		CoocMutAnalyzer an = new CoocMutAnalyzer(PortalDataset.BREAST_MUT);

		an.scoreReactionCoMutations();
		an.sortReacToScore();
		an.writeReactionsWithScore("/home/ozgun/Desktop/SortedReac.txt");

		Histogram hh = an.getReacScoreDistribution();

		Histogram his = new Histogram(1);
		int times = 100;
		Progress p = new Progress(times);
		for (int i = 0; i < times; i++)
		{
			p.tick();
			an.resetScores();
			an.shuffleGeneIdentitiesInReac();
			an.scoreReactionCoMutations();
			Histogram h = an.getReacScoreDistribution();
			his.add(h);
		}

		his.printTogether(hh, times);
	}

	public static void main(String[] args) throws IOException
	{
		PortalDataset dataset = PortalDataset.BREAST_MUT;
		CoocMutAnalyzer an = new CoocMutAnalyzer(dataset);

		an.prepareInteractingGenes();
		double[] scores = an.getUnexpectedDistribution(new HPRD());
//		double[] scores = an.getAllPairsDistribution();
		System.out.println("scores[0] = " + scores[0]);
		System.out.println("scores[1] = " + scores[1]);
		System.out.println("ratio = " + scores[0] / scores[1]);
		if (true) return;
		an.scoreInteractionCoMutations(new HPRD());
		ScoreUtil scTest = an.getInterScores();

		Histogram hh = an.getInterScoreHisto();

		Histogram his = new Histogram(1);
		int times = 10;
		Progress p = new Progress(times);
		ScoreUtil scRand = null;
		for (int i = 0; i < times; i++)
		{
			p.tick();
			an.resetScores();
			an.scoreInteractionCoMutations(new ShuffledHPRD());
			ScoreUtil sc = an.getInterScores();
			if (scRand == null) scRand = sc; else scRand.unite(sc);
			Histogram h = an.getInterScoreHisto();
			his.add(h);
		}

		double thr = scRand.getThresholdForFDR(scTest, 0.05);
		System.out.println("thr = " + thr);
		double fdr = scRand.getFDRForThr(scTest, thr);
		System.out.println("fdr = " + fdr);

		System.out.println("prev fdr = " + scRand.getFDRForThr(scTest, thr-1));

		his.printTogether(hh, times);

		an.scoreInteractionCoMutations(new HPRD());
		an.writeCoocInterAsSIF("C:\\Projects\\temp\\" + dataset.name() + ".sif", thr);
	}
}
