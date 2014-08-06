package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.model.RPPAData;
import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.network.SignedPC;
import org.cbio.causality.signednetwork.SignedType;
import org.cbio.causality.util.*;
import org.junit.Ignore;
import org.junit.Test;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class AnilsAnalysisCompareToCocktail
{
	public static final String DIR = "Anil-Data/";
	public static final double CHANGE_THR = 0.5;

	@Test
	@Ignore
	public void analyze() throws IOException
	{
		Map<String, RPPAData> probeMap = readABData();
		Map<String, Map<String, Set<RPPAData>>> ts = loadTreatments(probeMap);

		for (String treat : ts.keySet())
		{
			writeGraph(treat, ts.get(treat));
		}
	}

	public void writeGraph(String treat, Map<String, Set<RPPAData>> map) throws IOException
	{
		Set<String> sifLines = new HashSet<String>();
		Map<String, List<Relation>> relationMap = new HashMap<String, List<Relation>>();
		for (String time : map.keySet())
		{
			Set<RPPAData> set = map.get(time);
			selectChanged(set);
			List<Relation> rels = RPPANetworkMapper.map(set);
			RPPANetworkMapper.removeConflicting(rels);
			for (Relation rel : rels)
			{
				sifLines.add(rel.getEdgeData());
			}
			relationMap.put(time, rels);
		}

		BufferedWriter writer = new BufferedWriter(new FileWriter(DIR + treat + ".sif"));
		for (String line : sifLines) writer.write(line + "\n");
		writer.close();

		// prepare .formatseries file

		List<String> times = new ArrayList<String>(map.keySet());
		Collections.sort(times, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Double(o1.substring(0, o1.length() - 2)).compareTo(
					new Double(o2.substring(0, o2.length() - 2)));
			}
		});

		writer = new BufferedWriter(new FileWriter(DIR + treat + ".formatseries"));

		for (String time : times)
		{
			writer.write("group-name\t" + time + "\n");
			writer.write("node\tall-nodes\tcolor\t255 255 255\n");
			writer.write("node\tall-nodes\tbordercolor\t200 200 200\n");
			writer.write("node\tall-nodes\tborderwidth\t1\n");
			writer.write("node\tall-nodes\ttextcolor\t200 200 200\n");
			writer.write("edge\tall-edges\tcolor\t200 200 200\n");

			for (Relation rel : relationMap.get(time))
			{
				writer.write("edge\t" + rel.source + " " + rel.edgeType.getTag() + " " + rel.target +
					"\tcolor\t" + getEdgeColor(rel.edgeType) + "\n");
			}

			for (RPPAData data : map.get(time))
			{
				for (String gene : data.genes)
				{
					writer.write("node\t" + gene + "\ttextcolor\t0 0 0\n");
					if (data.isPhospho())
					{
						boolean unkwnEff = data.effect == null || data.effect == RPPAData.SiteEffect.COMPLEX;
						Color col = (unkwnEff && data.getChangeSign() > 0) ||
							(!unkwnEff && data.getActvityChangeSign() > 0) ?
							new Color(200, 100, 0) :
							new Color(50, 150, 200);

						writer.write("node\t" + gene + "\tbordercolor\t" +
							getColor(data.getChangeValue(), col) + "\n");

						if (!unkwnEff) writer.write("node\t" + gene + "\tborderwidth\t2\n");
					}
					else
					{
						Color col = data.getChangeSign() > 0 ? new Color(200, 100, 0) :
							new Color(50, 150, 200);

						writer.write("node\t" + gene + "\tcolor\t" +
							getColor(data.getChangeValue(), col) + "\n");
					}
				}
			}
		}
		writer.close();
	}

	private final static double MAX_VAL = 2;

	private String getColor(double val, Color maxCol)
	{
		val = Math.abs(val);
		double ratio = val / MAX_VAL;
		if (ratio > 1) ratio = 1;

		return (255 - (int) Math.round(ratio * (255 - maxCol.getRed()))) + " " +
			(255 - (int) Math.round(ratio * (255 - maxCol.getGreen()))) + " " +
			(255 - (int) Math.round(ratio * (255 - maxCol.getBlue())));
	}

	private String getEdgeColor(SignedType type)
	{
		switch (type)
		{
			case PHOSPHORYLATES:
			case UPREGULATES_EXPRESSION: return "0 100 0";
			case DEPHOSPHORYLATES:
			case DOWNREGULATES_EXPRESSION: return "100 0 0";
			default: return null;
		}
	}

	private void selectChanged(Set<RPPAData> set)
	{
		Iterator<RPPAData> iter = set.iterator();
		while (iter.hasNext())
		{
			RPPAData data = iter.next();
			if (Math.abs(data.getLog2MeanVal()) < CHANGE_THR) iter.remove();
		}
	}

	public Map<String, RPPAData> readABData() throws FileNotFoundException
	{
		Map<String, RPPAData> dataMap = new HashMap<String, RPPAData>();
		Scanner sc = new Scanner(new File(DIR + "ab_index_korkut.txt"));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			Probe p = new Probe(line);
			RPPAData data = new RPPAData(p.id, null, p.genes, p.sites);
			data.effect = p.activity == null ? null : p.activity ? RPPAData.SiteEffect.ACTIVATING :
				RPPAData.SiteEffect.INHIBITING;
			dataMap.put(p.id, data);
		}
		return dataMap;
	}

	private Map<String, Map<String, Set<RPPAData>>> loadTreatments(
		Map<String, RPPAData> probeMap) throws FileNotFoundException
	{
		Map<String, Map<String, Set<RPPAData>>> data = new HashMap<String, Map<String, Set<RPPAData>>>();

		Scanner sc = new Scanner(new File(DIR + "a2058_replicates.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] col = line.split("\t");
			String[] split = col[0].split("_");
			String treat = split[4];
			String time = split[5];

			String id = col[1];
			RPPAData rppa = (RPPAData) probeMap.get(id).clone();
			rppa.setChDet(new RPPAData.ChangeDetector()
			{
				@Override
				public int getChangeSign(RPPAData data)
				{
					return (int) Math.signum(getChangeValue(data));
				}

				@Override
				public double getChangeValue(RPPAData data)
				{
					return data.getLog2MeanVal();
				}
			});
			rppa.vals = new double[1][];
			rppa.vals[0] = new double[]{Double.parseDouble(col[5]), Double.parseDouble(col[6])};

			if (!data.containsKey(treat)) data.put(treat, new HashMap<String, Set<RPPAData>>());
			if (!data.get(treat).containsKey(time)) data.get(treat).put(time, new HashSet<RPPAData>());
			data.get(treat).get(time).add(rppa);
		}

		return data;
	}


	private class Probe
	{
		String id;
		List<String> genes;
		Map<String, List<String>> sites;
		boolean ph;
		Boolean activity;

		Probe(String line)
		{
			String[] split = line.split("\t");

			id = split[0];

			genes = new ArrayList<String>();
			Collections.addAll(genes, split[1].split("/"));

			ph = !split[2].equals("T");

			if (ph)
			{
				sites = new HashMap<String, List<String>>();
				List<String> list = new ArrayList<String>();
				for (String gene : genes)
				{
					sites.put(gene, list);
				}
				new HashSet<String>();
				split[2] = split[2].substring(split[2].indexOf("-") + 1, split[2].lastIndexOf(")"));
				for (String s : split[2].split(",/"))
				{
					String aa3 = s.substring(0, 3);
					String aa1 = aa3.equals("Tyr") ? "Y" : aa3.equals("Ser") ? "S" : aa3.equals("Thr") ? "T" : aa3.equals("Asp") ? "N" : "";

					if (aa1.isEmpty()) throw new RuntimeException(line);

					list.add(aa1 + s.substring(3));
				}

				boolean active = false;
				boolean inactive = false;

				for (String gene : genes)
				{
					for (String site : list)
					{
						Integer effect = PhosphoSitePlus.getEffect(gene, site);
						if (effect != null)
						{
							if (effect == 1) active = true;
							else if (effect == -1) inactive = true;
						}
					}
				}

				if (split[3].equals("a")) activity = true;
				else if (split[3].equals("i")) activity = false;

				Boolean pspAc = null;
				if (active && !inactive) pspAc = true;
				else if (!active && inactive) pspAc = false;

				if (ph && pspAc != null && !pspAc.equals(activity))
				{
					System.out.println("Mismatch in activity: " + line);
					System.out.println("pspAc = " + pspAc);
				}
			}
		}
	}

	@Test
	public void printData() throws FileNotFoundException
	{
		String[] content = new String[]{"901", "RB1"};

		Scanner sc = new Scanner(new File(DIR + "a2058_replicates.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			boolean select = true;
			for (String c : content)
			{
				if (!line.contains(c)) select = false;
			}
			if (select) System.out.println(line);
		}

	}


}
