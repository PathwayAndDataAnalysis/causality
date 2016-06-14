package org.cbio.causality.data.tcgafile;

import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Summary;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class MutationReader
{
	private Map<String, Map<String, String>> typeMap;
	private Map<String, Map<String, String>> valueMap;

	private Set<String> sampleSet;

	public static final String NO_DATA = "NO_DATA";

	public MutationReader(String filename) throws FileNotFoundException
	{
		this(filename, null);
	}

	public MutationReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this.typeMap = new LinkedHashMap<String, Map<String, String>>();
		this.valueMap = new LinkedHashMap<String, Map<String, String>>();
		this.sampleSet = new HashSet<String>();
		if (filename != null) load(filename, genes);
	}

	public void load(String filename, Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));

		int typeInd = -1;
		int sampleInd = -1;
		int protChInd = -1;

		Map<String, Map<String, Set<String>>> values = new HashMap<String, Map<String, Set<String>>>();

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("#")) continue;

			if (line.startsWith("Hugo_Symbol"))
			{
				String[] header = line.split("\t");
				typeInd = indexOf(header, "Variant_Classification");
				sampleInd = indexOf(header, "Tumor_Sample_Barcode");
				protChInd = indexOf(header, "Protein_Change");
				if (protChInd < 0) protChInd = indexOf(header, "amino_acid_change_WU");
				if (protChInd < 0) protChInd = indexOf(header, "AAChange");
				if (protChInd < 0) protChInd = indexOf(header, "amino_acid_change");

				if (protChInd < 0)
				{
					System.out.println("No protein change in file " + filename);
					return;
				}
				continue;
			}

			String id = line.substring(0, line.indexOf("\t"));

			if (genes != null && !genes.contains(id)) continue;

			String[] token = line.split("\t");

			String sample = token[sampleInd];
			sample = sample.substring(0, 15);

			sampleSet.add(sample);

			String type = token[typeInd];

			if (type.equals("Silent")) continue;

			String protCh = token.length <= protChInd ? "" : token[protChInd];
			if (!protCh.isEmpty() && protCh.startsWith("p.")) protCh = protCh.substring(2);
			else if (protCh.equals(".")) protCh = "";
			else if (!protCh.equals("NULL") && !protCh.isEmpty())
			{
//				System.out.println("protCh = " + protCh + " file=" + filename);
				continue;
			}

			if (!typeMap.containsKey(id)) typeMap.put(id, new HashMap<String, String>());
			typeMap.get(id).put(sample, type);

			if (!protCh.equals("NULL") && !protCh.isEmpty())
			{
				if (!values.containsKey(id)) values.put(id, new HashMap<String, Set<String>>());
				if (!values.get(id).containsKey(sample)) values.get(id).put(sample, new HashSet<String>());
				values.get(id).get(sample).add(protCh);
			}

			if (genes != null && genes.size() == typeMap.size()) break;
		}

		for (String gene : values.keySet())
		{
			if (!valueMap.containsKey(gene)) valueMap.put(gene, new HashMap<String, String>());
			for (String sample : values.get(gene).keySet())
			{
				if (valueMap.get(gene).containsKey(sample)) continue;

				StringBuilder s = new StringBuilder();
				for (String mut : values.get(gene).get(sample))
				{
					s.append(mut).append(" ");
				}
				valueMap.get(gene).put(sample, s.toString().trim());
			}
		}
	}

	public Set<String> getSamples()
	{
		return sampleSet;
	}

	public Set<String> getGenes()
	{
		return typeMap.keySet();
	}

	private int indexOf(String[] array, String val)
	{
		for (int i = 0; i < array.length; i++)
		{
			if (array[i].equals(val)) return i;
		}
		return -1;
	}

	/**
	 * All samples have to be in this dataset. This method does not support NO_DATA conditions.
	 */
	public boolean[] getGeneAlterationArray(String id, String[] samples)
	{
		if (typeMap.containsKey(id))
		{
			boolean[] b = new boolean[samples.length];
			Arrays.fill(b, false);
			for (int i = 0; i < samples.length; i++)
			{
				if (typeMap.get(id).containsKey(samples[i])) b[i] = true;
				else throw new IllegalArgumentException("Sample " + samples[i] + " does not have mutation data.");
			}
			return b;
		}
		return null;
	}

	public String[] getMutationValues(String id, String[] samples)
	{
		if (valueMap.containsKey(id))
		{
			String[] vals = new String[samples.length];
			for (int i = 0; i < samples.length; i++)
			{
				if (valueMap.get(id).containsKey(samples[i]))
				{
					vals[i] = valueMap.get(id).get(samples[i]);
				} else if (!sampleSet.contains(samples[i]))
				{
					vals[i] = NO_DATA;
				}
			}
			return vals;
		}
		return null;
	}

	public void writeAsAlterationMatrix(String outFile) throws IOException
	{
		List<String> samples = new ArrayList<String>(getSamples());
		Collections.sort(samples);
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (String gene : typeMap.keySet())
		{
			writer.write("\n" + gene);
			Map<String, String> map = typeMap.get(gene);

			for (String sample : samples)
			{
				writer.write("\t" + (map.containsKey(sample) ? "1" : "0"));
			}
		}

		writer.close();
	}

	private void printRecurrenceCounts()
	{
		int totalMut = 0;
		int delMut = 0;
		Map<String, Map<String, Integer>> cnt = new HashMap<String, Map<String, Integer>>();
		for (String gene : valueMap.keySet())
		{
			for (String sample : valueMap.get(gene).keySet())
			{
				for (String mut : valueMap.get(gene).get(sample).split(" "))
				{
					totalMut++;
					if (mut.contains("*") || mut.contains("fs")) delMut++;
					if (!cnt.containsKey(gene)) cnt.put(gene, new HashMap<String, Integer>());
					if (cnt.get(gene).containsKey(mut)) cnt.get(gene).put(mut, cnt.get(gene).get(mut) + 1);
					else cnt.get(gene).put(mut, 1);
				}
			}
		}

		System.out.println("Global ratio of deleterious mutations = " + (delMut / (double) totalMut));

		final Map<String, Integer> best = new HashMap<String, Integer>();
		for (String gene : cnt.keySet())
		{
			best.put(gene, Summary.max(cnt.get(gene).values()));
		}

		List<String> genes = new ArrayList<String>(cnt.keySet());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return best.get(o2).compareTo(best.get(o1));
			}
		});

		Map<String, Double> dRat = getRatiosOfDeleteriousMutations();

		Histogram h = new Histogram(0.1);
		h.setBorderAtZero(true);

		for (String gene : genes)
		{
			Integer c = best.get(gene);
			if (c == 2) break;
			System.out.println(c + "\t" + gene + "\t" + dRat.get(gene));
			h.count(dRat.get(gene));
		}
		h.print();
	}

	public Map<String, Integer> getHighestRecurrenceCounts()
	{
		Map<String, Map<String, Integer>> cnt = new HashMap<String, Map<String, Integer>>();
		for (String gene : valueMap.keySet())
		{
			for (String sample : valueMap.get(gene).keySet())
			{
				for (String mut : valueMap.get(gene).get(sample).split(" "))
				{
					if (!cnt.containsKey(gene)) cnt.put(gene, new HashMap<String, Integer>());
					if (cnt.get(gene).containsKey(mut)) cnt.get(gene).put(mut, cnt.get(gene).get(mut) + 1);
					else cnt.get(gene).put(mut, 1);
				}
			}
		}

		Map<String, Integer> highest = new HashMap<String, Integer>();
		for (String gene : cnt.keySet())
		{
			highest.put(gene, Summary.max(cnt.get(gene).values()));
		}
		return highest;
	}

	public Map<String, Double> getRatiosOfDeleteriousMutations()
	{
		Map<String, Double> rat = new HashMap<String, Double>();
		for (String gene : valueMap.keySet())
		{
			int total = 0;
			int del = 0;
			for (String sample : valueMap.get(gene).keySet())
			{
				String s = valueMap.get(gene).get(sample);
				total++;
				if (s.contains("*") || s.contains("fs")) del++;
			}
			double r = del / (double) total;
			rat.put(gene, r);
		}
		return rat;
	}

	public double getOverallDelMutRatio()
	{
		int total = 0;
		int del = 0;

		for (String gene : valueMap.keySet())
		{
			for (String sample : valueMap.get(gene).keySet())
			{
				String s = valueMap.get(gene).get(sample);
				total++;
				if (s.contains("*") || s.contains("fs")) del++;
			}
		}
		double r = del / (double) total;
		return r;
	}



	public Map<String, Integer> getMutatedSampleCounts()
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();
		for (String gene : valueMap.keySet())
		{
			cnt.put(gene, valueMap.get(gene).keySet().size());
		}
		return cnt;
	}

	public static void main(String[] args) throws IOException
	{
//		String dir = "/home/ozgun/Documents/TCGA/PanCan/";
//		MutationReader reader = new MutationReader(dir + "tcga_pancancer_082115.vep.filter_whitelisted.maf");
//		reader.writeAsAlterationMatrix(dir + "DataMatrix.txt");

		MutationReader reader = new MutationReader("/home/babur/Documents/TCGA/SKCM/mutation.maf");
//		reader.printRecurrenceCounts();
	}
}
