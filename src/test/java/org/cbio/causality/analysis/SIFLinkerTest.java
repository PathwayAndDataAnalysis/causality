package org.cbio.causality.analysis;

import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SIFLinkerTest
{
	@Test
	public void testLinking()
	{
		SIFLinker linker = new SIFLinker();
		linker.load("/home/ozgun/Projects/mutex/src/main/resources/network.txt");

		Set<String> from = new HashSet<String>(Arrays.asList("TP53", "PTEN", "CTNNB1", "PIK3CA", "ARID1A", "KRAS", "ARID5A"));

		List<String> set = linker.link(from, from, 0);

		for (String s : set)
		{
			System.out.println(s);
		}
	}

	@Test
	@Ignore
	public void prepareGunnarFile() throws IOException
	{
		String dir = "/home/ozgun/Desktop/gunnar/";

		Traverse trav = new Traverse();
		trav.load(dir + "SIF_uniprot.txt", Collections.<String>emptySet(),
			new HashSet<String>(Arrays.asList(
				"state-change", "degrades", "transactivate", "transinhibit")));

		String[] source = new String[]{
			dir + "gencodeV14.v7.pancan_subset.ensembleID.uniprot.curated.noTFnoSF.list",
			dir + "splicing_factors.ensembleID.uniprot.curated.noTF.list",
			dir + "transcription_factors.ensembleID.curated.uniprot.list"
		};

		for (String sourceName : source)
		{
			String outName = sourceName.substring(
				0, sourceName.lastIndexOf(".")) + "-WithUpstream.list";

			BufferedWriter writer = new BufferedWriter(new FileWriter(outName));

			Scanner sc = new Scanner(new File(sourceName));

			boolean firstLine = true;
			while (sc.hasNextLine())
			{
				String seed = sc.nextLine();

				if (!firstLine) writer.write("\n");
				else firstLine = false;

				writer.write(seed);

				Set<String> upstr = trav.goBFS(
					Collections.singleton(seed), Collections.singleton(seed), false);

				List<String> list = new ArrayList<String>(upstr);
				Collections.sort(list);

				for (String s : list)
				{
					writer.write("\t" + s);
				}
			}
			sc.close();
			writer.close();

		}
	}
}
