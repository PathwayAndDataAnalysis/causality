package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.model.RPPAData;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.PhosphoSitePlus;
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
public class AnilsAnalysisJQ1Round2
{
	public static final String DIR = "Anil-Data/JQ1/";
	public static final double FDR_CUTOFF = 0.05;

	@Test
	@Ignore
	public void analyze() throws IOException
	{
		Map<String, RPPAData> probeMap = readABData();
		loadTreatments(probeMap);
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
				List<String> headerList = new ArrayList<String>();
				List<Double> valuesList = new ArrayList<Double>();

				List<String> doses = new ArrayList<String>(map.get(celline).get(ab).keySet());
				Collections.sort(doses, new Comparator<String>()
				{
					@Override
					public int compare(String o1, String o2)
					{
						return new Double(o1).compareTo(new Double(o2));
					}
				});


				for (String dose : doses)
				{
					List<Double> vals = map.get(celline).get(ab).get(dose);
					int i = 0;
					for (Double val : vals)
					{
						headerList.add("dose-" + dose + "-rep-" + (++i));
						valuesList.add(val);
					}
				}

				RPPAData rppa = (RPPAData) probeMap.get(ab).clone();

				rppa.vals = new double[1][];
				rppa.vals[0] = new double[valuesList.size()];
				rppa.header = new String[1][rppa.vals[0].length];
				for (int i = 0; i < valuesList.size(); i++)
				{
					rppa.vals[0][i] = valuesList.get(i);
					rppa.header[0][i] = headerList.get(i);
				}
				data.get(celline).add(rppa);
			}
		}

		for (String celline : map.keySet())
		{
			RPPAData.write(data.get(celline), DIR + celline + "-transformed-data.txt");
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
			Collections.addAll(genes, split[2].split("\\|"));

			ph = !split[4].equals("T");

			if (ph)
			{
				sites = new HashMap<String, List<String>>();
				List<String> list = Arrays.asList(split[4].split("_"));
				for (String gene : genes)
				{
					sites.put(gene, list);
				}

				boolean active = false;
				boolean inactive = false;

				for (String gene : genes)
				{
					for (String site : sites.get(gene))
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
	@Ignore
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

	// Section: Difference of differences

	@Test
	@Ignore
	public void generateForDifferenceOfDifferences() throws IOException
	{
		List<RPPAData> igrov1 = RPPAFileReader.readAnnotation(DIR + "IGROV1-transformed-data.txt",
			"ID", "Symbol", "Site", "Effect");
		RPPAFileReader.addValues(igrov1, DIR + "IGROV1-transformed-data.txt", "ID",
			Arrays.asList("dose-0.01-rep-1\tdose-0.01-rep-2\tdose-0.01-rep-3\tdose-0.02-rep-1\tdose-0.02-rep-2\tdose-0.02-rep-3\tdose-0.04-rep-1\tdose-0.04-rep-2\tdose-0.04-rep-3\tdose-0.08-rep-1\tdose-0.08-rep-2\tdose-0.08-rep-3\tdose-0.16-rep-1\tdose-0.16-rep-2\tdose-0.16-rep-3".split("\t")),
			Arrays.asList("dose-0.31-rep-1\tdose-0.31-rep-2\tdose-0.31-rep-3\tdose-0.63-rep-1\tdose-0.63-rep-2\tdose-0.63-rep-3\tdose-1.25-rep-1\tdose-1.25-rep-2\tdose-1.25-rep-3\tdose-2.5-rep-1\tdose-2.5-rep-2\tdose-2.5-rep-3\tdose-5-rep-1\tdose-5-rep-2\tdose-5-rep-3".split("\t")));


		List<RPPAData> cov318 = RPPAFileReader.readAnnotation(DIR + "COV318-transformed-data.txt",
			"ID", "Symbol", "Site", "Effect");
		RPPAFileReader.addValues(cov318, DIR + "COV318-transformed-data.txt", "ID",
			Arrays.asList("dose-0.01-rep-1\tdose-0.01-rep-2\tdose-0.01-rep-3\tdose-0.02-rep-1\tdose-0.02-rep-2\tdose-0.02-rep-3\tdose-0.04-rep-1\tdose-0.04-rep-2\tdose-0.04-rep-3\tdose-0.08-rep-1\tdose-0.08-rep-2\tdose-0.08-rep-3\tdose-0.16-rep-1\tdose-0.16-rep-2\tdose-0.16-rep-3".split("\t")),
			Arrays.asList("dose-0.31-rep-1\tdose-0.31-rep-2\tdose-0.31-rep-3\tdose-0.63-rep-1\tdose-0.63-rep-2\tdose-0.63-rep-3\tdose-1.25-rep-1\tdose-1.25-rep-2\tdose-1.25-rep-3\tdose-2.5-rep-1\tdose-2.5-rep-2\tdose-2.5-rep-3\tdose-5-rep-1\tdose-5-rep-2\tdose-5-rep-3".split("\t")));

		Map<String, RPPAData> igrovMap = getMap(igrov1);

		double threshold = 0.05;
		ChDet det = new ChDet();
		det.setThreshold(threshold);

		for (RPPAData dCov : cov318)
		{
			RPPAData dIg = igrovMap.get(dCov.id);
			double[][] v = dCov.vals;
			dCov.vals = new double[4][];
			dCov.vals[0] = v[0];
			dCov.vals[1] = v[1];
			dCov.vals[2] = dIg.vals[0];
			dCov.vals[3] = dIg.vals[1];

			dCov.setChDet(det);
		}

		RPPANetworkMapper.writeGraph(cov318, -Math.log(threshold) / Math.log(2), DIR + "COV318-dif-of-dif.sif",
			RPPANetworkMapper.GraphType.COMPATIBLE_WITH_SITE_MATCH, null);
	}

	private Map<String, RPPAData> getMap(Collection<RPPAData> datas)
	{
		Map<String, RPPAData> map = new HashMap<String, RPPAData>();
		for (RPPAData data : datas)
		{
			map.put(data.id, data);
		}
		return map;
	}

	class ChDet extends RPPAData.ChangeAdapter
	{
		@Override
		public int getChangeSign(RPPAData data)
		{
			double pv = getPval(data);
			if (pv > threshold) return 0;
			double diff = getDiff(data);
			return diff > 0 ? 1 : -1;
		}

		@Override
		public double getChangeValue(RPPAData data)
		{
			double pv = getPval(data);
			double diff = getDiff(data);
			double val = -Math.log(pv) / Math.log(2);
			if (diff < 0) val *= -1;
			return val;
		}

		private double getOffset(RPPAData data)
		{
			return Summary.mean(data.vals[3]) - Summary.mean(data.vals[2]);
		}

		private double getPval(RPPAData data)
		{
			double offset = getOffset(data);
			double[] v = new double[data.vals[1].length];
			for (int i = 0; i < v.length; i++)
			{
				v[i] = data.vals[1][i] - offset;
			}

			return StudentsT.getPValOfMeanDifference(data.vals[0], v);
		}

		private double getDiff(RPPAData data)
		{
			double offset = getOffset(data);
			double[] v = new double[data.vals[1].length];
			for (int i = 0; i < v.length; i++)
			{
				v[i] = data.vals[1][i] - offset;
			}

			return Summary.mean(v) - Summary.mean(data.vals[0]);
		}

	}

	@Test
	@Ignore
	public void printNeighorhood()
	{
		Graph graph = PathwayCommons.getGraph(SIFEnum.CONTROLS_PHOSPHORYLATION_OF);
		graph.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
		Set<String> seed = new HashSet<String>();
		Collections.addAll(seed, ("YWHAE\n" +
			"EIF4EBP1\n" +
			"EIF4EBP1\n" +
			"TP53BP1\n" +
			"ACACA ACACB\n" +
			"ACACA\n" +
			"NCOA3\n" +
			"AKT1 AKT2 AKT3\n" +
			"AKT1 AKT2 AKT3\n" +
			"AKT1 AKT2 AKT3\n" +
			"CTNNA1\n" +
			"PRKAA1\n" +
			"PRKAA1\n" +
			"ANXA1\n" +
			"AR\n" +
			"BAD\n" +
			"BAK1\n" +
			"BAX\n" +
			"BCL2\n" +
			"BCL2L1\n" +
			"BCL2L1\n" +
			"BECN1\n" +
			"CTNNB1\n" +
			"BID\n" +
			"BCL2L11\n" +
			"JUN\n" +
			"KIT\n" +
			"MET\n" +
			"MET\n" +
			"MYC\n" +
			"RAF1\n" +
			"RAF1\n" +
			"CASP3\n" +
			"CASP7\n" +
			"CASP9\n" +
			"CAV1\n" +
			"PECAM1\n" +
			"CDC2\n" +
			"CHEK1\n" +
			"CHEK1\n" +
			"CHEK2\n" +
			"CHEK2\n" +
			"BIRC2\n" +
			"COL6A1\n" +
			"CLDN7\n" +
			"PTGS2\n" +
			"CCNB1\n" +
			"CCND1\n" +
			"CCNE1\n" +
			"PARK7\n" +
			"DVL3\n" +
			"CDH1\n" +
			"EEF2\n" +
			"EEF2K\n" +
			"EGFR\n" +
			"EGFR\n" +
			"EGFR\n" +
			"EGFR\n" +
			"EGFR\n" +
			"EIF4E\n" +
			"ESR1\n" +
			"ESR1\n" +
			"ERCC1\n" +
			"PTK2\n" +
			"FN1\n" +
			"FOXO3\n" +
			"FOXO3\n" +
			"GATA3\n" +
			"GSK3A GSK3B\n" +
			"GSK3A GSK3B\n" +
			"GSK3A GSK3B\n" +
			"ERBB2\n" +
			"ERBB2\n" +
			"ERBB3\n" +
			"ERBB3\n" +
			"IGF1R\n" +
			"IGFBP2\n" +
			"INPP4B\n" +
			"IRS1\n" +
			"MAPK8\n" +
			"MAPK9\n" +
			"KRAS\n" +
			"MAPK1 MAPK3\n" +
			"MAP2K1\n" +
			"MAP2K1\n" +
			"ERRFI1\n" +
			"MRE11A\n" +
			"MSH2\n" +
			"MSH6\n" +
			"MTOR\n" +
			"MTOR\n" +
			"CDH2\n" +
			"NFKB1\n" +
			"NF2\n" +
			"NOTCH1\n" +
			"NOTCH3\n" +
			"CDH3\n" +
			"CDKN1A\n" +
			"CDKN1B\n" +
			"CDKN1B\n" +
			"CDKN1B\n" +
			"MAPK14\n" +
			"MAPK14\n" +
			"TP53\n" +
			"RPS6KB1\n" +
			"RPS6KB1\n" +
			"RPS6KA1\n" +
			"PARP1\n" +
			"PXN\n" +
			"PCNA\n" +
			"PDK1\n" +
			"PIK3CA\n" +
			"PIK3R1\n" +
			"PRKCA\n" +
			"PRKCA\n" +
			"PGR\n" +
			"AKT1S1\n" +
			"PTCH1\n" +
			"PTEN\n" +
			"RAB11A\n" +
			"RAB25\n" +
			"RAD50\n" +
			"RAD51\n" +
			"RB1\n" +
			"RB1\n" +
			"RPS6\n" +
			"RPS6\n" +
			"DIABLO\n" +
			"SMAD1\n" +
			"SMAD3\n" +
			"SMAD4\n" +
			"SNAI2\n" +
			"SRC\n" +
			"SRC\n" +
			"SRC\n" +
			"STAT3\n" +
			"STAT5A\n" +
			"STMN1\n" +
			"SYK\n" +
			"MAPT\n" +
			"WWTR1\n" +
			"TFF1\n" +
			"C12ORF5\n" +
			"TSC2\n" +
			"VASP\n" +
			"KDR\n" +
			"XIAP\n" +
			"XRCC1\n" +
			"YAP1\n" +
			"YAP1\n" +
			"YBX1\n" +
			"YBX1\n").split("\\s+"));

		Set<String> neighbors = graph.getUpstream(seed);
		neighbors.addAll(seed);
		List<String> neigh = new ArrayList<String>(neighbors);
		Collections.sort(neigh);
		System.out.println("neigh = " + neigh);
		System.out.println("neigh.size() = " + neigh.size());
	}
}
