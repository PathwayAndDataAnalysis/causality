package org.cbio.causality.rppa;

import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.util.BaseDir;

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
	 * @param siteMatchStrict option to enforce matching a phosphorylation site in the network with
	 *                       the annotation of antibody
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @param customNetworkDirectory The directory that the network will be downloaded and SignedPC
	 *                               directory will be created in. Pass null to use default.
	 * @throws IOException
	 */
	public static void generateRPPAGraphs(String platformFile, String idColumn,
		String symbolsColumn, String sitesColumn, String effectColumn, String valuesFile,
		String valueColumn, double valueThreshold, String graphType, boolean siteMatchStrict,
		String outputFilePrefix, String customNetworkDirectory) throws IOException
	{
		if (customNetworkDirectory != null) BaseDir.setDir(customNetworkDirectory);

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
			siteMatchStrict ? RPPANetworkMapper.GraphType.CONFLICTING_WITH_SITE_MATCH : RPPANetworkMapper.GraphType.CONFLICTING :
			siteMatchStrict ? RPPANetworkMapper.GraphType.COMPATIBLE_WITH_SITE_MATCH : RPPANetworkMapper.GraphType.COMPATIBLE;

		// Generate output
		RPPANetworkMapper.writeGraph(datas, valueThreshold, outputFilePrefix + ".sif", type , null);
	}

	// Test in class. Bad practice. Tsk tsk tsk
	public static void main(String[] args) throws IOException
	{
		generateRPPAGraphs("/home/ozgun/Documents/JQ1/abdata-chibe.txt", "ID1", "Symbols", "Sites",
			"Effect", "/home/ozgun/Documents/JQ1/ovcar4_dif_drug_sig.txt", "change", 0.001,
			"compatible", true, "/home/ozgun/Temp/temp", null);
	}
}
