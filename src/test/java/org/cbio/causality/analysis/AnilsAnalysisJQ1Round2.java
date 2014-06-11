package org.cbio.causality.analysis;

import org.cbio.causality.model.RPPAData;
import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.signednetwork.SignedType;
import org.cbio.causality.util.CollectionUtil;
import org.cbio.causality.util.FDR;
import org.junit.Ignore;
import org.junit.Test;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class AnilsAnalysisJQ1Round2
{
	public static final String DIR = "Anil-Data/JQ1/";
	public static final double FDR_CUTOFF = 0.05;

	@Test
	@Ignore
	public void analyze() throws IOException
	{
		Map<String, RPPAData> probeMap = readABData();
		Map<String, Set<RPPAData>> ts = loadTreatments(probeMap);

		for (String celline : ts.keySet())
		{
			writeGraph(celline, ts.get(celline));
		}
	}

	public void writeGraph(String cellline, Set<RPPAData> datas) throws IOException
	{
		Set<String> sifLines = new HashSet<String>();
		selectChanged(datas);
		List<Relation> rels = RPPANetworkMapper.map(datas);
		RPPANetworkMapper.removeConflicting(rels);
		for (Relation rel : rels)
		{
			sifLines.add(rel.getEdgeData());
		}

		BufferedWriter writer = new BufferedWriter(new FileWriter(DIR + cellline + ".sif"));
		for (String line : sifLines) writer.write(line + "\n");
		writer.close();

		// prepare .formatseries file

		writer = new BufferedWriter(new FileWriter(DIR + cellline + ".format"));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t200 200 200\n");
		writer.write("node\tall-nodes\tborderwidth\t1\n");
		writer.write("node\tall-nodes\ttextcolor\t200 200 200\n");
		writer.write("edge\tall-edges\tcolor\t200 200 200\n");

		for (Relation rel : rels)
		{
			writer.write("edge\t" + rel.source + " " + rel.edgeType.getTag() + " " + rel.target +
				"\tcolor\t" + getEdgeColor(rel.edgeType) + "\n");
		}

		for (RPPAData data : datas)
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

		Map<String, Set<String>> gene2rppa = new HashMap<String, Set<String>>();
		for (RPPAData data : datas)
		{
			for (String gene : data.genes)
			{
				if (!gene2rppa.containsKey(gene)) gene2rppa.put(gene, new HashSet<String>());
				gene2rppa.get(gene).add(data.id);
			}
		}

		for (String gene : gene2rppa.keySet())
		{
			writer.write("node\t" + gene + "\ttooltip\t" +
				CollectionUtil.merge(gene2rppa.get(gene), "\\n") + "\n");
		}
		writer.close();
	}

	private final static double MAX_VAL = 15;

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
		Map<String, Double> pvals = new HashMap<String, Double>();
		for (RPPAData data : set)
		{
			pvals.put(data.id, data.getTTestPval());
		}
		List<String> select = FDR.selectBH(pvals, FDR_CUTOFF);
		Iterator<RPPAData> iter = set.iterator();
		while (iter.hasNext())
		{
			RPPAData data = iter.next();
			if (!select.contains(data.id)) iter.remove();
		}
	}

	public Map<String, RPPAData> readABData() throws FileNotFoundException
	{
		Map<String, RPPAData> dataMap = new HashMap<String, RPPAData>();
		Scanner sc = new Scanner(new File(DIR + "abdata.txt"));
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

	private Map<String, Set<RPPAData>> loadTreatments(Map<String, RPPAData> probeMap) throws FileNotFoundException
	{
		Map<String, Set<RPPAData>> data = new HashMap<String, Set<RPPAData>>();

		Scanner sc = new Scanner(new File(DIR + "data-merged.txt"));

		String header = sc.nextLine();
		String[] abname = header.substring(header.indexOf("\t1") + 1).split("\t");

		Map<String, Map<String, Map<String, List<Double>>>> map =
			new HashMap<String, Map<String, Map<String, List<Double>>>>();

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] col = line.split("\t");

			String[] split = col[0].split("_");
			String cellline = split[1];
			String dose = split[3];

			if (dose.equals("10")) continue;

			for (int i = 2; i < col.length; i++)
			{
				double val = Double.parseDouble(col[i]);

				if (!map.containsKey(cellline)) map.put(cellline, new HashMap<String, Map<String, List<Double>>>());
				if (!map.get(cellline).containsKey(abname[i - 2])) map.get(cellline).put(abname[i - 2], new HashMap<String, List<Double>>());
				if (!map.get(cellline).get(abname[i - 2]).containsKey(dose)) map.get(cellline).get(abname[i - 2]).put(dose, new ArrayList<Double>());
				map.get(cellline).get(abname[i - 2]).get(dose).add(val);
			}
		}

		for (String celline : map.keySet())
		{
			if (!data.containsKey(celline)) data.put(celline, new HashSet<RPPAData>());

			for (String ab : map.get(celline).keySet())
			{
				List<Double> lowVals = new ArrayList<Double>();
				List<Double> highVals = new ArrayList<Double>();

				for (String dose : map.get(celline).get(ab).keySet())
				{
					List<Double> vals = map.get(celline).get(ab).get(dose);
					if (Double.parseDouble(dose) > 0.2) highVals.addAll(vals);
					else lowVals.addAll(vals);
				}

				RPPAData rppa = (RPPAData) probeMap.get(ab).clone();
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
						return data.getSignificanceBasedVal();
					}
				});
				rppa.vals0 = new double[lowVals.size()];
				for (int i = 0; i < lowVals.size(); i++)
				{
					rppa.vals0[i] = lowVals.get(i);
				}
				rppa.vals1 = new double[highVals.size()];
				for (int i = 0; i < highVals.size(); i++)
				{
					rppa.vals1[i] = highVals.get(i);
				}
				data.get(celline).add(rppa);
			}
		}

		return data;
	}


	private class Probe
	{
		String id;
		Set<String> genes;
		Set<String> sites;
		boolean ph;
		Boolean activity;

		Probe(String line)
		{
			String[] split = line.split("\t");

			id = split[0];

			genes = new HashSet<String>();
			Collections.addAll(genes, split[2].split("\\|"));

			ph = !split[4].equals("T");

			if (ph)
			{
				sites = new HashSet<String>(Arrays.asList(split[4].split("_")));

				boolean active = false;
				boolean inactive = false;

				for (String gene : genes)
				{
					for (String site : sites)
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

					activity = pspAc;
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
