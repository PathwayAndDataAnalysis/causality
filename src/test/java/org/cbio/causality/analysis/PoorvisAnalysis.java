package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.rppa.RPPAData;
import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.network.SignedPC;
import org.cbio.causality.signednetwork.SignedType;
import org.cbio.causality.util.*;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
@Ignore
public class PoorvisAnalysis
{
	public static final String DIR = "Poorvi-Data/";
	public static final boolean PVAL_BY_PERMUTATION = false;
	public static final int ITERATION = 10000;

	@Test
	@Ignore
	public void analyze() throws IOException
	{
		Map<String, Treatment> ts = loadTreatments();

		System.out.println("probes = " + ts.values().iterator().next().valsMap.size());

		Map<String, Probe> probeMap = readGeneData();

		int cntp = 0;
		int cnta = 0;
		for (Probe probe : probeMap.values())
		{
			if (probe.ph)
			{
				if (probe.activity == null) cnta++;
				else cntp++;
			}
		}
		System.out.println("unkwn phospho effect = " + cnta);
		System.out.println("known phospho effect = " + cntp);

		for (String drug : ts.keySet())
		{
			Treatment tr = ts.get(drug);
//			if (!tr.name.equals("RO-31-7549")) continue;

			Map<String, Double> pvals = new HashMap<String, Double>();
			Map<Probe, int[]> windows = new HashMap<Probe, int[]>();
			for (String probeID : tr.valsMap.keySet())
			{
				Probe probe = probeMap.get(probeID);
				windows.put(probe, tr.searchBestWindow(probeID));
				pvals.put(probeID, tr.getPvalOfWindow(probeID, windows.get(probe)));
			}
			List<String> selected = FDR.select(pvals, null, 0.05);
			System.out.println(tr.name + "\t" + selected.size());

			List<Probe> diffProbes = new ArrayList<Probe>(selected.size());
			for (String probeID : selected)
			{
				assert probeMap.containsKey(probeID) : probeID;
				diffProbes.add(probeMap.get(probeID));
			}

			for (int i = 0; i < tr.times.size(); i++)
			{
				System.out.println(i + "\t" +
					countNonZero(tr.getProbeStatus(diffProbes, windows, i)));
			}

			for (String probeID : selected)
			{
				Probe probe = probeMap.get(probeID);
				int[] win = windows.get(probe);
				System.out.print((win[1] - win[0] + 1) + "\t" + Arrays.toString(win));
				System.out.print(tr.getDirection(probeID, win) ? "\tup" : "\tdw");
				System.out.print(probe.ph ? (probe.activity == null ? "\tukw" : probe.activity ?
					"\tact" : "\tinh") : "\t");
				System.out.println("\t" + probeID + " " + probe.genes + "\t" +
					FormatUtil.roundToSignificantDigits(pvals.get(probeID), 2));
			}

			Map<Probe, Integer> status = tr.getProbeStatus(diffProbes, windows);

			Graph[] graph = new Graph[SignedType.values().length];
			Map<String[], int[]>[] edge2window = new Map[graph.length];
			int i = 0;
			for (SignedType type : SignedType.values())
			{
				Map<String[], int[]> e2w = new HashMap<String[], int[]>();
				graph[i] = getGraph(type, diffProbes, status, e2w, windows);
				edge2window[i++] = e2w;
			}
			String siffile = DIR + drug + ".sif";
			write(graph, siffile);
			String formatfile = DIR + drug + ".formatseries";
			writeSeriesFormat(tr, graph, edge2window, formatfile, diffProbes, windows);
			checkSIFSanity(siffile, formatfile);
		}
	}

	private void write(Graph[] graph, String filename) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (Graph g : graph)
		{
			g.write(writer);
		}

		writer.close();
	}

	private int countNonZero(Map<Probe, Integer> status)
	{
		int cnt = 0;
		for (Probe probe : status.keySet())
		{
			if (status.get(probe) != 0) cnt++;
		}
		return cnt;
	}

	public Map<String, Probe> readGeneData() throws FileNotFoundException
	{
		Map<String, Probe> map = new HashMap<String, Probe>();
		Scanner sc = new Scanner(new File(DIR + "AntibodyData.txt"));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			Probe p = new Probe(line);
			map.put(p.id, p);
		}
		return map;
	}

	private Map<String, Treatment> loadTreatments() throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(DIR + "TS676_SingleDrugTimeSeries_2014March03.txt"));
		String header = sc.nextLine();

		String[] split = header.split("\t");
		List<String> temp = Arrays.asList(split).subList(11, split.length);
		List<String> genes = new ArrayList<String>(temp.size());
		for (String s : temp)
		{
			if (s.startsWith("X") && Character.isDigit(s.charAt(1))) s = s.substring(1);
			genes.add(s);
		}

		Map<String, Treatment> treatments = new HashMap<String, Treatment>();
		List<String> lines = new ArrayList<String>();
		String current = null;

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String name = line.split("\t")[3];

			if (current != null && (!current.equals(name) || !sc.hasNextLine()))
			{
				if (!sc.hasNextLine()) lines.add(line);
				Treatment tr = new Treatment(lines, genes);
				treatments.put(tr.name, tr);
				lines.clear();
			}

			lines.add(line);
			current = name;
		}
		return treatments;
	}

	private Graph getGraph(SIFType type, List<Probe> probes,
		Map<Probe, Integer> status, Map<String[], int[]> edge2window, Map<Probe, int[]> windows)
	{
		Graph graph = SignedPC.getGraph(type);
		Graph g2 = new Graph(type.getTag(), graph.getEdgeType());

		for (Probe probe1 : probes)
		{
			for (Probe probe2 : probes)
			{
				if (probe1 == probe2) continue;

				int[] ov = getWindowOverlap(probe1, probe2, windows);
				if (ov == null || windows.get(probe2)[0] < ov[0]) continue;

				List<String[]> edges = getEdges(probe1, probe2, graph);

				if (edges.isEmpty()) continue;
				boolean posRel = status.get(probe1) * status.get(probe2) *
					probe1.getAssumedActivity() == 1;

				if (type == SignedType.PHOSPHORYLATES)
				{
					if (!posRel) continue;
					if (!probe2.ph) continue;
				}
				else if (type == SignedType.DEPHOSPHORYLATES)
				{
					if (posRel) continue;
					if (!probe2.ph) continue;
				}
				else if (type == SignedType.UPREGULATES_EXPRESSION)
				{
					if (!posRel) continue;
					if (probe2.ph) continue;
				}
				else if (type == SignedType.DOWNREGULATES_EXPRESSION)
				{
					if (posRel) continue;
					if (probe2.ph) continue;
				}

				for (String[] edge : edges)
				{
					g2.putRelation(edge[0], edge[1], true);
					edge2window.put(edge, ov);
				}
			}
		}
		return g2;
	}

	private int[] getWindowOverlap(Probe p1, Probe p2, Map<Probe, int[]> windows)
	{
		int[] w1 = windows.get(p1);
		int[] w2 = windows.get(p2);
		int[] ov = new int[2];
		ov[0] = Math.max(w1[0], w2[0]);
		ov[1] = Math.min(w1[1], w2[1]);

		if (ov[0] <= ov[1]) return ov;
		return null;
	}

	private List<String[]> getEdges(Probe p1, Probe p2, Graph graph)
	{
		List<String[]> list = new ArrayList<String[]>();

		for (String gene1 : p1.genes)
		{
			Set<String> down = graph.getDownstream(gene1);
			for (String gene2 : p2.genes)
			{
				if (down.contains(gene2))
				{
					list.add(new String[]{gene1, gene2});
				}
			}
		}
		return list;
	}

	private void writeSeriesFormat(Treatment tr, Graph[] graph, Map<String[], int[]>[] e2w,
		String filename, List<Probe> probes, Map<Probe, int[]> windows)
	{
		try
		{
			Set<String> genes = new HashSet<String>();
			for (Graph g : graph) genes.addAll(g.getSymbols());

			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

			for (int i = 0; i < tr.times.size(); i++)
			{
				writer.write("group-name\t" + tr.times.get(i) + " hr\n");

				writer.write("node\tall-nodes\tcolor\t255 255 255\n");
				writer.write("node\tall-nodes\tbordercolor\t200 200 200\n");
				writer.write("node\tall-nodes\tborderwidth\t1\n");
				writer.write("node\tall-nodes\ttextcolor\t200 200 200\n");

				writer.write("edge\tall-edges\tcolor\t200 200 200\n");

				Set<Probe> onPrb = new HashSet<Probe>();
				Map<Probe, Integer> status = tr.getProbeStatus(probes, windows, i);
				for (Probe p1 : probes)
				{
					if (status.get(p1) != 0)
					{
						onPrb.add(p1);
						String color = status.get(p1) == 1 ? "200 255 200" : "255 255 200";
						String bordercolor = status.get(p1) * p1.getAssumedActivity() == 1 ?
							"250 150 100" : "100 200 250";
						for (String gene : p1.genes)
						{
							if (genes.contains(gene))
							{
								writer.write("node\t" + gene + "\ttextcolor\t0 0 0\n");

								if (p1.ph)
								{
									writer.write("node\t" + gene + "\tborderwidth\t2\n");
									writer.write("node\t" + gene + "\tbordercolor\t" + bordercolor +
										"\n");
								}
								else
								{
									writer.write("node\t" + gene + "\tcolor\t" + color + "\n");
								}
							}
						}
					}
				}

				for (Probe p1 : onPrb)
				{
					for (String target : p1.genes)
					{
						for (int j = 0; j < graph.length; j++)
						{
							Set<String> up = graph[j].getUpstream(target);

							for (Probe p2 : onPrb)
							{
								if (p1 == p2) continue;

								for (String source : p2.genes)
								{
									if (up.contains(source) && inWindow(source, target, i, e2w[j]))
									{
										writer.write("edge\t" + source + " " +
											graph[j].getEdgeType() + " " + target + "\tcolor\t" +
											getEdgeActiveColor(graph[j].getEdgeType()) + "\n");
									}
								}
							}

						}
					}
				}
			}

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private boolean inWindow(String source, String target, int win, Map<String[], int[]> e2w)
	{
		for (String[] edge : e2w.keySet())
		{
			if (edge[0].equals(source) && edge[1].equals(target))
			{
				int[] w = e2w.get(edge);

				if (win >= w[0] && win <= w[1]) return true;
			}
		}
		return false;
	}

	private static Map<String, String> edge2color;

	private String getEdgeActiveColor(String type)
	{
		if (edge2color == null)
		{
			edge2color = new HashMap<String, String>();
			edge2color.put(SignedType.PHOSPHORYLATES.getTag(), "0 100 0");
			edge2color.put(SignedType.DEPHOSPHORYLATES.getTag(), "100 0 0");
			edge2color.put(SignedType.UPREGULATES_EXPRESSION.getTag(), "0 80 20");
			edge2color.put(SignedType.DOWNREGULATES_EXPRESSION.getTag(), "80 20 0");
		}
		return edge2color.get(type);
	}

	private void checkSIFSanity(String siffile, String formatfile) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(siffile));
		Set<String> edges = new HashSet<String>();
		while(sc.hasNextLine())
		{
			edges.add(sc.nextLine().replaceAll("\t", " "));
		}

		sc = new Scanner(new File(formatfile));
		Set<String> found = new HashSet<String>();
		while(sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("edge\t"))
			{
				String[] split = line.split("\t");
				found.add(split[1]);
			}
		}

		System.out.println("edges = " + edges.size());
		System.out.println("used = " + found.size());

		found.removeAll(edges);
		for (String s : found)
		{
			System.out.println(s);
		}
		System.out.println();
	}

	class Treatment
	{
		String name;
		List<Integer> times;
		int[] timeArray;
		Map<String, double[]> valsMap;

		Treatment(List<String> lines, List<String> probes)
		{
			times = new ArrayList<Integer>();
			timeArray = new int[lines.size()];
			valsMap = new HashMap<String, double[]>();

			for (String probe : probes)
			{
				valsMap.put(probe, new double[lines.size()]);
			}

			int i = 0;
			for (String line : lines)
			{
				String[] split = line.split("\t");

				if (name == null) name = split[3];
				timeArray[i] = Integer.parseInt(split[5]);

				if (!times.contains(timeArray[i])) times.add(timeArray[i]);

				for (int j = 0; j < probes.size(); j++)
				{
					valsMap.get(probes.get(j))[i] = Double.parseDouble(split[j + 11]);
				}

				i++;
			}
		}

		public double getPvalOfWindow(String gene, int... win)
		{
			boolean[][] m = getSubsetMarkers(putTimesInSet(win));
			double[][] v = getValues(gene, m);
			return PVAL_BY_PERMUTATION ?
				StudentsT.getPValOfMeanDifferenceBySimulation(v[0], v[1], ITERATION, 10) :
				StudentsT.getPValOfMeanDifference(v[0], v[1]);
		}

		public boolean getDirection(String probe, int[] win)
		{
			boolean[][] m = getSubsetMarkers(putTimesInSet(win));
			double[][] v = getValues(probe, m);
			return Summary.calcChangeOfMean(v[0], v[1]) > 0;
		}

		private Set<Integer> putTimesInSet(int[] win)
		{
			Set<Integer> timeSet = new HashSet<Integer>();
			for (int k = win[0]; k <= win[1]; k++)
			{
				timeSet.add(times.get(k));
			}
			return timeSet;
		}

		public int[] searchBestWindow(String probe)
		{
//			if (true) return new int[]{0, 2};

			double bestPval = 1;
			int[] win = new int[2];
			for (int width = times.size() - 1; width >= 1; width--)
			{
				for (int i = 1; i < times.size() - width + 1; i++)
				{
					int j = i + width - 1;

					double p = getPvalOfWindow(probe, i, j);

					if (p < bestPval)
					{
						bestPval = p;
						win[0] = i;
						win[1] = j;
					}
				}
			}
			return win;
		}

		private boolean[][] getSubsetMarkers(Set<Integer> timeSet)
		{
			boolean[][] m = new boolean[2][timeArray.length];
			int min = Summary.min(timeSet);

			for (int i = 0; i < timeArray.length; i++)
			{
				m[0][i] = timeArray[i] < min;
				m[1][i] = timeSet.contains(timeArray[i]);
			}
			return m;
		}

		private double[][] getValues(String probe, boolean[][] m)
		{
			double[][] v = new double[2][];
			v[0] = getValues(probe, m[0]);
			v[1] = getValues(probe, m[1]);
			return v;
		}

		private double[] getValues(String probe, boolean[] m)
		{
			double[] v = new double[ArrayUtil.countValue(m, true)];

			int j = 0;
			for (int i = 0; i < timeArray.length; i++)
			{
				if (m[i]) v[j++] = valsMap.get(probe)[i];
			}
			return v;
		}

		public Map<Probe, Integer> getProbeStatus(List<Probe> probes,
			Map<Probe, int[]> windows, int timeIndex)
		{
			Map<Probe, Integer> map = new HashMap<Probe, Integer>();

			for (Probe probe : probes)
			{
				int[] win = windows.get(probe);

				if (win[0] <= timeIndex && win[1] >= timeIndex)
				{
					map.put(probe, getDirection(probe.id, win) ? 1 : -1);
				}
				else map.put(probe, 0);
			}
			return map;
		}

		public Map<Probe, Integer> getProbeStatus(List<Probe> probes,
			Map<Probe, int[]> windows)
		{
			Map<Probe, Integer> map = new HashMap<Probe, Integer>();

			for (Probe probe : probes)
			{
				int[] win = windows.get(probe);
				map.put(probe, getDirection(probe.id, win) ? 1 : -1);
			}
			return map;
		}

	}

	class Probe
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
			id = id.replaceAll("-", ".").replaceAll(" ", ".");

			genes = new HashSet<String>();

			for (String term : split[2].split("\\."))
			{
				for (String gene : term.split("_"))
				{
					gene = HGNC.getSymbol(gene);
					if (gene != null) genes.add(gene);
				}
			}

			ph = split[3].equals("1");

			if (ph)
			{
				sites = new HashSet<String>();
				for (String s : id.replaceAll("_", " ").replaceAll("\\.", " ").split(" "))
				{
					if (s.startsWith("p")) sites.add(s.substring(1));
				}

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

				if (!active && !inactive)
				{
					for (String gene : genes)
					{
						for (String site : sites)
						{
							Integer effect = PhosphoSitePlus.getClosestEffect(gene, site);
							if (effect != null)
							{
								if (effect == 1) active = true;
								else if (effect == -1) inactive = true;
							}
						}
					}
				}

				if (active && !inactive) activity = true;
				else if (!active && inactive) activity = false;
			}
		}

		public int getAssumedActivity()
		{
			if (activity == null) return 1;
			else if (!activity) return -1;
			else return 1;
		}

		@Override
		public String toString()
		{
			return id;
		}
	}

	// Section: Data conversion

	@Test
	@Ignore
	public void convert() throws FileNotFoundException
	{
		Map<String, RPPAData> dataMap = readABFile(DIR + "AntibodyData.txt");
		Map<String, Set<RPPAData>> map = readTreatments(DIR +
			"TS676_SingleDrugTimeSeries_2014July31.txt", dataMap);

		augmentDMSO(map);

		for (String treatment : map.keySet())
		{
			RPPAData.write(map.get(treatment), DIR + "formatted/" + treatment + ".txt");
		}
	}

	Map<String, RPPAData> readABFile(String filename) throws FileNotFoundException
	{
		Map<String, RPPAData> map = new HashMap<String, RPPAData>();

		Scanner sc = new Scanner(new File(filename));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");

			String id = token[0].replaceAll(" ", ".").replaceAll("-", ".");
			
			List<String> genes = new ArrayList<String>();

			Map<String, List<String>> sitesMap = null;

			for (String gene : token[2].split("\\."))
			{
				if (gene.contains("_"))
				{
					String sites = gene.substring(gene.indexOf("_"));
					gene = gene.substring(0, gene.indexOf("_"));

					for (String site : sites.split("_"))
					{
						if (site.startsWith("p")) site = site.substring(1);
						if (!site.isEmpty())
						{
							if (sitesMap == null) sitesMap = new HashMap<String, List<String>>();
							if (!sitesMap.containsKey(gene)) sitesMap.put(gene, new ArrayList<String>());
							sitesMap.get(gene).add(site);
						}
					}
				}
				genes.add(gene);
			}
			
			RPPAData data = new RPPAData(id, null, genes, sitesMap);

			map.put(id, data);
		}
		return map;
	}

	Map<String, Set<RPPAData>> readTreatments(String filename, Map<String, RPPAData> dataMap) throws FileNotFoundException
	{
		String[] ids;
		Map<String, Set<RPPAData>> map = new HashMap<String, Set<RPPAData>>();

		String treatment = null;

		Scanner sc = new Scanner(new File(filename));
		String first = sc.nextLine();
		first = first.substring(first.indexOf("X"));
		ids = first.split("\t");
		for (int i = 0; i < ids.length; i++)
		{
			if (ids[i].startsWith("X")) ids[i] = ids[i].substring(1);
			else break;
		}

		Map<String, RPPAData> copy = copyMap(dataMap);
		Map<String, List<Double>> idToVals = new HashMap<String, List<Double>>();
		List<String> header = new ArrayList<String>();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");

			String treat = token[1];
			if (treatment != null && !treat.equals(treatment))
			{
				HashSet<RPPAData> set = new HashSet<RPPAData>(copy.values());
				associateValues(set, idToVals, header);
				map.put(treatment, set);
				copy = copyMap(dataMap);
				idToVals = new HashMap<String, List<Double>>();
				treatment = treat;
				header = new ArrayList<String>();
			}
			else if (treatment == null) treatment = treat;

			String rowID = token[2] + "_" + token[3];
			header.add(rowID);

			for (int i = 0; i < token.length - 4; i++)
			{
				double val = Double.parseDouble(token[i + 4]);
				if (!idToVals.containsKey(ids[i])) idToVals.put(ids[i], new ArrayList<Double>());
				idToVals.get(ids[i]).add(val);
			}
		}

		HashSet<RPPAData> set = new HashSet<RPPAData>(copy.values());
		associateValues(set, idToVals, header);
		map.put(treatment, set);

		return map;
	}

	private Map<String, RPPAData> copyMap(Map<String, RPPAData> map)
	{
		Map<String, RPPAData> copy = new HashMap<String, RPPAData>();
		for (String id : map.keySet())
		{
			copy.put(id, (RPPAData) map.get(id).clone());
		}
		return copy;
	}

	private void associateValues(Set<RPPAData> datas, Map<String, List<Double>> valMap,
		List<String> header)
	{
		int size = header.size();
		assert size == valMap.values().iterator().next().size();

		String[] h = new String[header.size()];
		for (int i = 0; i < size; i++)
		{
			h[i] = header.get(i);
		}

		for (RPPAData data : datas)
		{
			data.vals = new double[1][];
			data.vals[0] = new double[size];
			data.header = new String[][]{h};

			for (int i = 0; i < size; i++)
			{
				data.vals[0][i] = valMap.get(data.id).get(i);
			}
		}
	}

	private void augmentDMSO(Map<String, Set<RPPAData>> dataMap)
	{
		Set<RPPAData> dmso = dataMap.get("DMSO");
		Set<RPPAData> dmsoH = dataMap.get("DMSOh");

		for (String treatment : dataMap.keySet())
		{
			if (treatment.equals("DMSO") || treatment.equals("DMSOh")) continue;

			if (treatment.equals("Temozolomide")) augmentDMSO(dataMap.get(treatment), dmsoH);
			else augmentDMSO(dataMap.get(treatment), dmso);
		}
		dataMap.remove("DMSO");
		dataMap.remove("DMSOh");
	}

	private void augmentDMSO(Set<RPPAData> treatData, Set<RPPAData> dmsoData)
	{
		Map<String, RPPAData> tMap = putInMap(treatData);
		Map<String, RPPAData> dMap = putInMap(dmsoData);

		String[] header = null;

		for (String id : tMap.keySet())
		{
			RPPAData tData = tMap.get(id);
			RPPAData dData = dMap.get(id);

			if (header == null)
			{
				header = new String[tData.header[0].length + dData.header[0].length];
				for (int i = 0; i < dData.header[0].length; i++)
				{
					header[i] = "dmso_" + dData.header[0][i];
				}
				for (int i = 0; i < tData.header[0].length; i++)
				{
					header[i + dData.header[0].length] = "treat_" + tData.header[0][i];
				}
			}

			tData.header[0] = header;
			double[] vals = new double[header.length];
			System.arraycopy(dData.vals[0], 0, vals, 0, dData.vals[0].length);
			System.arraycopy(tData.vals[0], 0, vals, dData.vals[0].length, tData.vals[0].length);
			tData.vals[0] = vals;
		}
	}

	private Map<String, RPPAData> putInMap(Set<RPPAData> datas)
	{
		Map<String, RPPAData> map = new HashMap<String, RPPAData>();
		for (RPPAData data : datas)
		{
			map.put(data.id, data);
		}
		return map;
	}
}
