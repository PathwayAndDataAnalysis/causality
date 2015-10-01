package org.cbio.causality.data.tcgafile;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class MutationReader
{
	private String filename;

	private Map<String, Map<String, String>> data;

	public MutationReader(String filename) throws FileNotFoundException
	{
		this(filename, null);
	}
	public MutationReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this.filename = filename;
		this.data = new LinkedHashMap<String, Map<String, String>>();
		load(genes);
	}

	private void load(Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		String line = sc.nextLine();
		while (line.startsWith("#")) line = sc.nextLine();

		String[] header = line.split("\t");

//		int idInd = indexOf(header, "Hugo_Symbol");
		int typeInd = indexOf(header, "Variant_Classification");
		int sampleInd = indexOf(header, "Tumor_Sample_Barcode");

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			if (line.startsWith("Hugo_Symbol")) continue;
			String id = line.substring(0, line.indexOf("\t"));

			if (genes != null && !genes.contains(id)) continue;

			String[] token = line.split("\t");
			String type = token[typeInd];

			if (type.equals("Silent")) continue;

			String sample = token[sampleInd];
			sample = sample.substring(0, 15);

			if (!data.containsKey(id)) data.put(id, new HashMap<String, String>());
			data.get(id).put(sample, type);

			if (genes != null && genes.size() == data.size()) break;
		}
	}

	public Set<String> getSamples()
	{
		Set<String> samples = new HashSet<String>();
		for (Map<String, String> map : data.values())
		{
			samples.addAll(map.keySet());
		}
		return samples;
	}

	public Set<String> getGenes()
	{
		return data.keySet();
	}

	private int indexOf(String[] array, String val)
	{
		for (int i = 0; i < array.length; i++)
		{
			if (array[i].equals(val)) return i;
		}
		return -1;
	}

	public boolean[] getGeneAlterationArray(String id, String[] samples)
	{
		boolean[] b = new boolean[samples.length];
		Arrays.fill(b, false);
		if (data.containsKey(id))
		{
			for (int i = 0; i < samples.length; i++)
			{
				if (data.get(id).containsKey(samples[i])) b[i] = true;
			}
		}
		return b;
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

		for (String gene : data.keySet())
		{
			writer.write("\n" + gene);
			Map<String, String> map = data.get(gene);

			for (String sample : samples)
			{
				writer.write("\t" + (map.containsKey(sample) ? "1" : "0"));
			}
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgun/Documents/TCGA/PanCan/";
		MutationReader reader = new MutationReader(dir + "tcga_pancancer_082115.vep.filter_whitelisted.maf");
		reader.writeAsAlterationMatrix(dir + "DataMatrix.txt");
	}
}
