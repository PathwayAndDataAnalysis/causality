package org.cbio.causality.data.tcgafile;

import org.cbio.causality.rppa.RPPAData;
import org.cbio.causality.util.CollectionUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RPPAReader
{
	private String filename;

	private static Map<String, Set<String>> additionalAnnotation;

	private Map<String, Map<String, Double>> data;

	private Map<String, Map<String, Set<String>>> symbolToIDs;

	private Map<String, RPPAData> idToData;

	public RPPAReader(String filename) throws FileNotFoundException
	{
		this(filename, null);
	}

	public RPPAReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this.filename = filename;
		this.data = new HashMap<String, Map<String, Double>>();
		this.symbolToIDs = new HashMap<String, Map<String, Set<String>>>();
		this.idToData = new HashMap<String, RPPAData>();
		load(genes);
	}

	private void load(Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		String line = sc.nextLine();
		while (line.startsWith("#")) line = sc.nextLine();

		String[] header = line.split("\t");

		int ss = 0;
		while (!header[ss].startsWith("TCGA")) ss++;

		for (int i = ss; i < header.length; i++)
		{
			header[i] = header[i].substring(0, 15);
		}

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String id = line.substring(0, line.indexOf("\t"));

			Set<String> pGenes = parseID(id);
			id = id.substring(id.indexOf("|") + 1);


			Map<String, String> map = mapSymbols(pGenes);

			if (genes != null && !CollectionUtil.intersects(map.keySet(), genes)) continue;

			List<String> geneList = new ArrayList<String>(map.keySet());
			Collections.sort(geneList);

			Map<String, List<String>> sitesMap = new HashMap<String, List<String>>();
			for (String g : geneList)
			{
				String pGene = map.get(g);
				if (!pGene.contains("_"))
				{
					sitesMap = null;
					break;
				}
				pGene = pGene.substring(pGene.indexOf("_") + 1);
				List<String> sites = new ArrayList<String>();
				Collections.addAll(sites, pGene.split("_"));
				sitesMap.put(g, sites);
			}

			RPPAData d = new RPPAData(id, null, geneList, sitesMap);

			idToData.put(id, d);

			for (String sym : map.keySet())
			{
				if (!symbolToIDs.containsKey(sym)) symbolToIDs.put(sym, new HashMap<String, Set<String>>());
				if (!symbolToIDs.get(sym).containsKey(map.get(sym)))
					symbolToIDs.get(sym).put(map.get(sym), new HashSet<String>());
				symbolToIDs.get(sym).get(map.get(sym)).add(id);
			}

			String[] token = line.split("\t");

			for (int i = ss; i < header.length; i++)
			{
				double val = token[i].equals("NA") ? Double.NaN : Double.parseDouble(token[i]);

				if (!data.containsKey(id)) data.put(id, new HashMap<String, Double>());
				data.get(id).put(header[i], val);
			}
		}
	}

	private Map<String, String> mapSymbols(Set<String> pGenes)
	{
		Map<String, String> map = new TreeMap<String, String>();
		for (String pGene : pGenes)
		{
			String symbol = pGene.contains("_") ? pGene.substring(0, pGene.indexOf("_")) : pGene;
			map.put(symbol, pGene);
		}
		return map;
	}

	public Set<String> getSamples()
	{
		Set<String> samples = new HashSet<String>();
		for (Map<String, Double> map : data.values())
		{
			samples.addAll(map.keySet());
		}
		return samples;
	}

	public Set<String> getGenes()
	{
		return symbolToIDs.keySet();
	}

	public double[] getValues(String id, String[] samples)
	{
		double[] b = new double[samples.length];
		if (data.containsKey(id))
		{
			for (int i = 0; i < samples.length; i++)
			{
				if (data.get(id).containsKey(samples[i]))
				{
					b[i] = data.get(id).get(samples[i]);
				}
				else
				{
					b[i] = Double.NaN;
				}
			}
		}
		return b;
	}

	private Set<String> parseID(String id)
	{
		String symbols = id.substring(0, id.indexOf("|")).trim();
		id = id.substring(id.indexOf("|") + 1).trim();

		if (additionalAnnotation.containsKey(id)) return additionalAnnotation.get(id);

		List<String> pSites = new ArrayList<String>();

		if (id.contains("_p"))
		{
			String s = id.substring(id.indexOf("_p") + 2);

			String[] token = s.split("_");

			for (int i = 0; i < token.length; i++)
			{
				if (token[i].startsWith("p")) token[i] = token[i].substring(1);
				pSites.add(token[i]);
			}
		}

		Set<String> genes = new HashSet<String>(Arrays.asList(symbols.split(" ")));

		if (!pSites.isEmpty() && genes.size() > 1)
		{
			System.out.println("Warning: There is an unhandled RPPA entry: " + symbols + "|" + id);
		}

		Set<String> pGenes = new HashSet<String>();

		for (String gene : genes)
		{
			for (String site : pSites)
			{
				gene += "_" + site;
			}
			pGenes.add(gene);
		}
		return pGenes;
	}

	public Set<RPPAData> getAssociatedData(String symbol, String[] samples)
	{
		Map<String, Set<String>> map = symbolToIDs.get(symbol);
		if (map == null) return Collections.emptySet();

		Set<RPPAData> stubs = new HashSet<RPPAData>();

		for (String pGene : map.keySet())
		{
			for (String id : map.get(pGene))
			{
				stubs.add(idToData.get(id));
			}
		}

		Set<RPPAData> set = new HashSet<RPPAData>();

		for (RPPAData stub : stubs)
		{
			RPPAData d = (RPPAData) stub.clone();
			d.vals = new double[][]{getValues(d.id, samples)};
			set.add(d);
		}

		return set;
	}

	static
	{
		additionalAnnotation = new HashMap<String, Set<String>>();
		Scanner sc = new Scanner(RPPAReader.class.getResourceAsStream("rppa-annotation.txt"));
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");

			additionalAnnotation.put(token[0], new HashSet<String>(Arrays.asList(token).subList(1, token.length)));
		}
		sc.close();
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		RPPAReader reader = new RPPAReader("/home/babur/Documents/Temp/BRCA/rppa.txt");
		List<String> sList = new ArrayList<String>(reader.getSamples());
		Collections.sort(sList);
		String[] samples = sList.toArray(new String[sList.size()]);
		for (RPPAData data : reader.getAssociatedData("AKT1", samples))
		{
			System.out.println(data.toString());
		}
	}
}