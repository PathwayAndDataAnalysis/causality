package org.cbio.causality.data.tcgafile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ExpressionReader
{
	private String filename;

	private Map<String, Map<String, Double>> data;

	private double LOG2 = Math.log(2);

	public ExpressionReader(String filename) throws FileNotFoundException
	{
		this (filename, null);
	}

	public ExpressionReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this.filename = filename;
		this.data = new HashMap<String, Map<String, Double>>();
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

		// skip second line
		sc.nextLine();

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String id = line.substring(0, line.indexOf("|"));

			if (genes != null && !genes.contains(id)) continue;
			if (id.equals("?")) continue;

			String[] token = line.split("\t");

			for (int i = ss; i < header.length; i++)
			{
				Double val = Double.parseDouble(token[i]);

				if (!data.containsKey(id)) data.put(id, new HashMap<String, Double>());
				data.get(id).put(header[i], val);
			}

			if (genes != null && data.size() == genes.size()) break;
		}
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
		return data.keySet();
	}

	public double[] getGeneAlterationArray(String id, String[] samples)
	{
		double[] d = new double[samples.length];

		if (data.containsKey(id))
		{
			for (int i = 0; i < samples.length; i++)
			{
				if (data.get(id).containsKey(samples[i])) d[i] = Math.log1p(data.get(id).get(samples[i])) / LOG2;
			}
		}
		return d;
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		ExpressionReader reader = new ExpressionReader("/home/ozgun/Documents/TCGA/UVM/UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt");
		System.out.println(reader.getSamples().size());
	}
}
