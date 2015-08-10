package org.cbio.causality.data.tcgafile;

import org.cbio.causality.util.Download;
import org.cbio.causality.util.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class BroadDownloader
{
	private static final String BROAD_URL_PREFIX = "http://gdac.broadinstitute.org/runs/";
	private static final String BROAD_DATA_URL_PREFIX = BROAD_URL_PREFIX + "stddata__";
	private static final String BROAD_ANALYSIS_URL_PREFIX = BROAD_URL_PREFIX + "analyses__";


	private static final String MUT_ARCH_PART = "/gdac.broadinstitute.org_?.Mutation_Packager_Calls.Level_3.";
	private static final String EXP_ARCH_PART = "/gdac.broadinstitute.org_?.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.";
	private static final String CNA_ARCH_PART = "/gdac.broadinstitute.org_?-TP.CopyNumber_Gistic2.Level_4.";

	//http://gdac.broadinstitute.org/runs/stddata__2015_04_02/data/GBM/201504022015040200.0.0.tar.gz

	public static void downloadAll(String date, String dir) throws IOException
	{
		List<String> codes = getStudyCodes(date);
//		codes.retainAll(Arrays.asList("FPPP"));

		for (String code : codes)
		{
			if (!download(date, dir, code))
			{
				FileUtil.delete(new File (dir + File.separator + code));
			}
		}
	}

	public static boolean download(String date, String dir, String code) throws IOException
	{
		return
			downloadCopyNumber(date, dir, code) &&
			downloadExpression(date, dir, code);// &&
//			downloadMutations(date, dir, code);
	}

	public static boolean downloadCopyNumber(String date, String dir, String code) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		new File(directory).mkdirs();
		String d = date.replaceAll("_", "");

		String url = BROAD_ANALYSIS_URL_PREFIX + date + "/data/" + code + "/" + d + CNA_ARCH_PART.replace("?",code) + d + "00.0.0.tar.gz";
		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "all_thresholded.by_genes", directory + "copynumber.txt"))
			{
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	public static boolean downloadExpression(String date, String dir, String code) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		new File(directory).mkdirs();
		String d = date.replaceAll("_", "");

		String url = BROAD_DATA_URL_PREFIX + date + "/data/" + code + "/" + d + EXP_ARCH_PART.replace("?",code) + d + "00.0.0.tar.gz";
		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractEntryContainingNameInTARGZFile(tempFile, "data.txt", directory + "expression.txt"))
			{
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	public static boolean downloadMutations(String date, String dir, String code) throws IOException
	{
		String directory = dir + File.separator + code + File.separator;
		new File(directory).mkdirs();
		String d = date.replaceAll("_", "");

		String url = BROAD_DATA_URL_PREFIX + date + "/data/" + code + "/" + d + MUT_ARCH_PART.replace("?",code) + d + "00.0.0.tar.gz";
		String tempFile = directory + "temp.tar.gz";
		if (Download.downloadAsIs(url, tempFile))
		{
			if (FileUtil.extractAllEntriesContainingNameInTARGZFile(tempFile, "maf.txt", directory))
			{
				uniteMutations(directory);
				return new File(tempFile).delete();
			}
		}
		return false;
	}

	private static void uniteMutations(String dir) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "mutation.maf"));
		File d = new File(dir);
		for (File file : d.listFiles())
		{
			if (file.getName().endsWith(".maf.txt"))
			{
				Scanner sc = new Scanner(file);
				while (sc.hasNextLine())
				{
					writer.write(sc.nextLine() + "\n");
				}
				file.delete();
			}
		}
		writer.close();
	}

	public static List<String> getStudyCodes(String date) throws IOException
	{
		List<String> codes = new ArrayList<String>();
		Scanner sc = new Scanner(new URL(BROAD_DATA_URL_PREFIX + date + "/ingested_data.tsv").openStream());
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("Totals")) break;

			codes.add(line.substring(0, line.indexOf("\t")));
		}
		return codes;
	}

	public static void main(String[] args) throws IOException
	{
//		List<String> codes = getStudyCodes("2015_06_01");
//		System.out.println(codes);

		download("2015_04_02", "/home/ozgun/Documents/TCGA/broad", "SARC");

//		downloadAll("2015_04_02", "/home/ozgun/Documents/TCGA/broad");
	}
}
