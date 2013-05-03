package org.cbio.causality.data.portal;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.cbio.causality.idmapping.EntrezGene;
import org.cbio.causality.idmapping.HGNC;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * @author Arman Aksoy
 * @author Ozgun Babur
 */

public class CBioPortalManager
{
	private static Log log = LogFactory.getLog(CBioPortalManager.class);

	private static final String CACHE_DIR = "portal-cache";

	private static String portalURL = "http://www.cbioportal.org/public-portal/webservice.do?";
	protected final static String COMMAND = "cmd=";
	protected final static String DELIMITER = "\t";

	public CBioPortalManager()
	{
	}

	public static String getPortalURL()
	{
		return CBioPortalManager.portalURL;
	}

	public static void setPortalURL(String portalURL)
	{
		CBioPortalManager.portalURL = portalURL;
	}

	protected String[] downloadDataForGene(String symbol, GeneticProfile geneticProfile, CaseList caseList)
		throws IOException
	{
		String geneid = EntrezGene.getID(symbol);
		String url = "getProfileData&case_set_id=" + caseList.getId() + "&"
			+ "genetic_profile_id=" + geneticProfile.getId() + "&"
			+ "gene_list=" + geneid;

		List<String[]> results = parseURL(url, false);

		if (results.size() < 2)
		{
			log.warn("Cannot get data for " + symbol);
			return null;
		}

		// DEBUG CODE
		String[] cases = results.get(0);
		for (int i = 0; i < caseList.getCases().length; i++)
		{
			assert cases[i+2].equals(caseList.getCases()[i]);
		}
		// END OF DEBUG CODE


		String[] data = results.get(1);

		assert data.length > 2;

		String[] result = new String[data.length - 2];
		System.arraycopy(data, 2, result, 0, result.length);
		return result;
	}

	public List<CaseList> getCaseListsForStudy(CancerStudy study) throws IOException
	{
		List<CaseList> caseLists = new ArrayList<CaseList>();

		String url = "getCaseLists&cancer_study_id=" + study.getStudyId();
		for (String[] results : parseURL(url))
		{
			assert results.length == 5;
			String[] cases = results[4].split(" ");
			assert cases.length > 0;

			CaseList caseList = new CaseList(results[0], results[1], cases);
			caseLists.add(caseList);
		}

		return caseLists;
	}

	private List<String[]> parseURL(String urlPostFix) throws IOException
	{
		return parseURL(urlPostFix, true);
	}

	private List<String[]> parseURL(String urlPostFix, boolean skipHeader) throws IOException
	{
		List<String[]> list = new ArrayList<String[]>();

		String urlStr = portalURL + COMMAND + urlPostFix;
		URL url = new URL(urlStr);
		URLConnection urlConnection = url.openConnection();
		Scanner scanner = new Scanner(urlConnection.getInputStream());

		int lineNum = 0;
		while (scanner.hasNextLine())
		{
			String line = scanner.nextLine();
			lineNum++;

			if (line.startsWith("#") || line.length() == 0 || (skipHeader && lineNum == 2))
				continue;

			list.add(line.split(DELIMITER));
		}

		return list;
	}

	public List<GeneticProfile> getGeneticProfilesForStudy(CancerStudy study) throws IOException
	{
		List<GeneticProfile> geneticProfiles = new ArrayList<GeneticProfile>();

		String url = "getGeneticProfiles" + "&cancer_study_id=" + study.getStudyId();
		for (String[] results : parseURL(url))
		{
			assert results.length == 6;
			GeneticProfile geneticProfile = new GeneticProfile(results[0], results[1], results[2], results[4]);
			geneticProfiles.add(geneticProfile);
		}

		assert !geneticProfiles.isEmpty();
		return geneticProfiles;
	}

	public List<CancerStudy> getCancerStudies() throws IOException
	{
		List<CancerStudy> studies = new ArrayList<CancerStudy>();
		String urlStr = "getCancerStudies";
		for (String[] result : parseURL(urlStr))
		{
			assert result.length == 3;
			CancerStudy cancerStudy = new CancerStudy(result[0], result[1], result[2]);
			studies.add(cancerStudy);
		}
		return studies;
	}

	public void cacheData(String[] data, String symbol, GeneticProfile geneticProfile,
		CaseList caseList)
	{
		String dir = CACHE_DIR + File.separator + geneticProfile.getId() + File.separator +
			caseList.getId();

		File f = new File(dir);
		if (!f.exists()) f.mkdirs();

		String url = dir  + File.separator + symbol;
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(url));

			for (int i = 0; i < data.length; i++)
			{
				writer.write(data[i]);
				if (i < data.length - 1) writer.write(DELIMITER);
			}

			writer.close();
		}
		catch (IOException e)
		{
			log.error("Cannot cache data for " + symbol, e);
		}
	}

	public String[] readDataInCache(String symbol, GeneticProfile geneticProfile, CaseList caseList)
	{
		String url = CACHE_DIR + File.separator + geneticProfile.getId() + File.separator +
			caseList.getId() + File.separator + symbol;

		if (new File(url).exists())
		{
			try
			{
				BufferedReader reader = new BufferedReader(new FileReader(url));
				String line = reader.readLine();
				reader.close();
				return line.split(DELIMITER);
			}
			catch (IOException e)
			{
				log.error("Cannot read an existing file", e);
			}
		}
		return null;
	}


	public String[] getDataForGene(String symbol, GeneticProfile geneticProfile, CaseList caseList)
	{
		try
		{
			String[] data = readDataInCache(symbol, geneticProfile, caseList);
			if (data != null) return data;

			data = downloadDataForGene(symbol, geneticProfile, caseList);
			if (data != null)
			{
				cacheData(data, symbol, geneticProfile, caseList);
			}
			return data;
		}
		catch (IOException e)
		{
			log.error("Cannot access to the data of " + symbol, e);
			return null;
		}
	}
}
