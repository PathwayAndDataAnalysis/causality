package org.cbio.causality.data.portal;

import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.Download;
import org.cbio.causality.util.FileUtil;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class BroadAccessor
{
	private static final String BROAD_DIR = "broad-data/";
	private static final String CACHED_STUDIES_FILE = "studies.txt";
	private static String cacheDir;
	private static String broadDataURL = "http://gdac.broadinstitute.org/runs/analyses__2013_05_23/";
	private static List<String> studyCodes;
	private static final String MUTSIG_ANALYSIS_SUBSTR = "MutSigNozzleReportMerged.Level_4";
	private static final String GISTIC_ANALYSIS_SUBSTR = "Gistic2.Level_4";

	public static void setCacheDir(String dir)
	{
		cacheDir = dir;
	}

	public static void setBroadDataURL(String url)
	{
		broadDataURL = url;
		studyCodes = null;
	}

	public static String getBroadDataURL()
	{
		return broadDataURL;
	}

	public static List<String> getStudyCodes()
	{
		if (studyCodes == null)
		{
			studyCodes = readStudiesFromCache();

			if (studyCodes == null)
			{
				studyCodes = new ArrayList<String>(30);
				try
				{
					URL url = new URL(broadDataURL + "ingested_data.tsv");

					URLConnection con = url.openConnection();

					BufferedReader reader = new BufferedReader(
						new InputStreamReader(con.getInputStream()));

					for (String line = reader.readLine(); line != null; line = reader.readLine())
					{
						if (line.isEmpty() || line.startsWith("#")
							|| line.startsWith("Tumor") || line.startsWith("Totals")) continue;

						String study = line.substring(0, line.indexOf("\t"));
						if (bothMutsigAndGisticAvailable(study)) studyCodes.add(study);
					}
					reader.close();

					// Keep only the ones that are available in cBioPortal

					Set<String> available = new HashSet<String>();
					CBioPortalAccessor acc = new CBioPortalAccessor();
					for (CancerStudy cancerStudy : acc.getCancerStudies())
					{
						if (cancerStudy.getStudyId().endsWith("_tcga"))
						{
							available.add(cancerStudy.getStudyId().substring(0,
								cancerStudy.getStudyId().indexOf("_")).toUpperCase());
						}
					}
					studyCodes.retainAll(available);
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}

				if (!studyCodes.isEmpty())
				{
					cacheStudies(studyCodes);
				}
			}
		}
		return studyCodes;
	}

	private static List<String> readStudiesFromCache()
	{
		try
		{
			if (!new File(getStudiesFileName()).exists()) return null;

			List<String> studies = new ArrayList<String>();
			BufferedReader reader = new BufferedReader(new FileReader(getStudiesFileName()));

			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				if (!line.isEmpty()) studies.add(line);
			}

			reader.close();
			return studies;
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	private static void cacheStudies(List<String> studies)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(
				new FileWriter(getStudiesFileName()));

			for (String study : studies)
			{
				writer.write(study + "\n");
			}

			writer.close();

		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private static String getStudiesFileName()
	{
		return getBroadCacheDir() + CACHED_STUDIES_FILE;
	}

	private static String getBroadDateString()
	{
		String s = broadDataURL;
		s = s.substring(s.indexOf("__") + 2, s.lastIndexOf("/"));
		s = s.replaceAll("_", "");
		return s;
	}

	private static String getBroadDataURL(String study)
	{
		return broadDataURL + "data/" + study + "/" + getBroadDateString() + "/";
	}

	private static String getBroadCacheDir()
	{
		if (cacheDir == null)
		{
			String s = CBioPortalManager.getCacheDir() + BROAD_DIR;
			File f = new File(s);
			if (!f.exists()) f.mkdirs();
			return s;
		}
		return cacheDir;
	}

	private static List<String> getBroadAnalysisFileNames(String study)
	{
		List<String> list = new ArrayList<String>(30);
		try
		{
			URL url = new URL(getBroadDataURL(study));

			URLConnection con = url.openConnection();
			BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()));
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				String start = "<li><a href=\"";
				if (line.startsWith(start))
				{
					String file = line.substring(start.length(), line.indexOf("\">"));
					list.add(file);
				}
			}
			reader.close();
		}
		catch (IOException e)
		{
			System.out.println(e);
		}
		return list;
	}

	private static String getGisticFileName(List<String> list)
	{
		for (String s : list)
		{
			if (s.contains(GISTIC_ANALYSIS_SUBSTR)) return s;
		}
		return null;
	}

	private static String getMutsigFileName(List<String> list)
	{
		for (String s : list)
		{
			if (s.contains(MUTSIG_ANALYSIS_SUBSTR)) return s;
		}
		return null;
	}

	private static boolean bothMutsigAndGisticAvailable(String study)
	{
		List<String> analysisFiles = getBroadAnalysisFileNames(study);
		return !analysisFiles.isEmpty() &&
			getGisticFileName(analysisFiles) != null &&
			getMutsigFileName(analysisFiles) != null;
	}

	private static String getCachedGisticFileName(String study)
	{
		return getBroadCacheDir() + study + "_gistic.tar.gz";
	}

	private static String getCachedMutsigFileName(String study)
	{
		return getBroadCacheDir() + study + "_mutsig.tar.gz";
	}

	private static boolean downloadGistic(String study, List<String> analysisFileNames)
	{
		String s = getGisticFileName(analysisFileNames);
		return s != null &&
			Download.downloadAsIs(getBroadDataURL(study) + s, getCachedGisticFileName(study));
	}

	private static boolean downloadMutsig(String study, List<String> analysisFileNames)
	{
		String s = getMutsigFileName(analysisFileNames);
		return s != null &&
			Download.downloadAsIs(getBroadDataURL(study) + s, getCachedMutsigFileName(study));
	}

	public static Set<String> getMutsigGenes(String study, double qvalThr)
	{
		Set<String> genes = new HashSet<String>();
		String file = getCachedMutsigFileName(study);
		if (!new File(file).exists())
		{
			downloadMutsig(study, getBroadAnalysisFileNames(study));
		}
		if (new File(file).exists())
		{
			genes.addAll(readGenesFromMutsig(file, qvalThr));
		}

		return genes;
	}

	public static Set<String> getGisticGenes(String study, double qvalThr)
	{
		Set<String> genes = new HashSet<String>();
		List<Set<String>> sets = getGisticGeneSets(study, qvalThr);
		for (Set<String> set : sets)
		{
			genes.addAll(set);
		}
		return genes;
	}

	public static List<Set<String>> getGisticGeneSets(String study, double qvalThr)
	{
		String file = getCachedGisticFileName(study);
		if (!new File(file).exists())
		{
			downloadGistic(study, getBroadAnalysisFileNames(study));
		}
		if (new File(file).exists())
		{
			return readGenesFromGistic(file, qvalThr);
		}
		return null;
	}

	private static List<Set<String>> readGenesFromGistic(String filename, double qvalThr)
	{
		List<Set<String>> list = new ArrayList<Set<String>>();
		String s = FileUtil.readEntryContainingNameInTARGZFile(filename, "amp_genes");
		readGisticData(list, s, qvalThr);

		s = FileUtil.readEntryContainingNameInTARGZFile(filename, "del_genes");
		readGisticData(list, s, qvalThr);
		return list;
	}

	private static Set<String> readGenesFromMutsig(String filename, double qvalThr)
	{
		Set<String> set = new HashSet<String>();
		String s = FileUtil.readEntryContainingNameInTARGZFile(filename, ".sig_genes");
		if (s == null) s = FileUtil.readEntryContainingNameInTARGZFile(filename, "sig_genes");

		for (String line : s.split("\n"))
		{
			if (line.startsWith("rank")) continue;

			String[] token = line.split("\t");

			double qval = Double.parseDouble(token[token.length - 1]);

			if (qval < qvalThr)
			{
				String symbol = HGNC.getSymbol(token[1]);
				if (symbol != null) set.add(symbol);
			}
		}

		System.out.println("mutsig set = " + set.size());
		return set;
	}

	private static void readGisticData(List<Set<String>> list, String s, double qvalThr)
	{
		String[] line = s.split("\n");

		String[] qvalStr = line[2].split("\t");
		double[] qvals = new double[qvalStr.length - 1];
		for (int i = 0; i < qvals.length; i++)
		{
			qvals[i] = Double.parseDouble(qvalStr[i + 1]);
		}

		Set<String>[] set = new Set[qvals.length];
		for (int i = 0; i < set.length; i++)
		{
			set[i] = new HashSet<String>();
		}

		for (int i = 4; i < line.length; i++)
		{
			String[] gene = line[i].split("\t");

			for (int j = 1; j < gene.length; j++)
			{
				if (qvals[j - 1] >= qvalThr) continue;

				String symbol = HGNC.getSymbol(gene[j]);
				if (symbol != null) set[j - 1].add(symbol);
			}
		}
		for (Set<String> mem : set)
		{
			if (!mem.isEmpty()) list.add(mem);
		}
	}
}
