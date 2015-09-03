package org.cbio.causality.rppa;

import org.cbio.causality.network.PhosphoSitePlus;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * This class reads the RPPA platform and data files, and generates a ChiBE SIF graph.
 *
 * @author Ozgun Babur
 */
public class RPPAFrontFace
{
	/**
	 * Reads the RPPA platform and data files, and generates a ChiBE SIF graph.
	 *
	 * @param platformFile Name of the antibody reference file
	 * @param idColumn Column name of IDs
	 * @param symbolsColumn Column name of gene symbols
	 * @param sitesColumn Column name for phosphorylation sites
	 * @param effectColumn Column name for effect of the site on activity
	 * @param valuesFile Name of the measurements file
	 * @param valueColumn Name of the values column in the measurements file
	 * @param valueThreshold The value threshold to be considered as significant
	 * @param graphType Either "compatible" or "conflicting"
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @throws IOException
	 */
	public static void generateRPPAGraphs(String platformFile, String idColumn,
		String symbolsColumn, String sitesColumn, String effectColumn, String valuesFile,
		String valueColumn, double valueThreshold, String graphType, String outputFilePrefix)
		throws IOException
	{
		// Read platform file
		List<RPPAData> datas = RPPAFileReader.readAnnotation(platformFile, idColumn, symbolsColumn,
			sitesColumn, effectColumn);

		// Read values
		List<String> vals0 = Arrays.asList(valueColumn);
		RPPAFileReader.addValues(datas, valuesFile, idColumn, vals0, null, 0D);

		// Prepare change detector class
		RPPAData.ChangeAdapter chDet = new RPPAData.ChangeAdapter(){};
		chDet.setThreshold(valueThreshold);
		for (RPPAData data : datas) if (!data.isActivity()) data.setChDet(chDet);

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.fillInMissingEffect(datas);

		// Set the graph type
		RPPANetworkMapper.GraphType type = graphType.toLowerCase().startsWith("conflict") ?
			RPPANetworkMapper.GraphType.CONFLICTING_WITH_SITE_MATCH :
			RPPANetworkMapper.GraphType.COMPATIBLE_WITH_SITE_MATCH;

		// Generate output
		RPPANetworkMapper.writeGraph(datas, valueThreshold, outputFilePrefix + ".sif", type , null);
	}

	// Test in class. Bad practice. Tsk tsk tsk
	public static void main(String[] args) throws IOException
	{
		generateRPPAGraphs("/home/ozgun/Documents/JQ1/abdata-chibe.txt", "ID1", "Symbols", "Sites",
			"Effect", "/home/ozgun/Documents/JQ1/ovcar4_dif_drug_sig.txt", "change", 0.001,
			"compatible", "/home/ozgun/Temp/temp");
	}
}
