package org.cbio.causality;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFSearcher;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SIFGenerator
{
	public static void main(String[] args) throws IOException
	{
		generate("/home/ozgun/Projects/biopax-pattern/All-Data.owl",
//		generate("/home/ozgun/Desktop/temp.owl",
			"SIF.txt");
	}

	public static void generate(String inputModelFile, String outputFile) throws IOException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream(inputModelFile));

		List<SIFInteraction> sifs = new ArrayList<SIFInteraction>(
			generate(model));

		Collections.sort(sifs);

		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

		for (SIFInteraction sif : sifs)
		{
			writer.write(sif + "\n");
		}

		writer.close();
	}

	public static Set<SIFInteraction> generate(Model model)
	{
		SIFSearcher searcher = new SIFSearcher(SIFType.CONTROLS_STATE_CHANGE,
			SIFType.CONTROLS_EXPRESSION, SIFType.CONTROLS_DEGRADATION);

		searcher.setMassDataMode(true);

		return searcher.searchSIF(model);
	}
}
