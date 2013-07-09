package org.cbio.causality.analysis;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Named;
import org.biopax.paxtools.model.level3.SmallMolecule;
import org.biopax.paxtools.model.level3.SmallMoleculeReference;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SIFLinkerTest
{
	@Test
	@Ignore
	public void testSmallMoelcules() throws FileNotFoundException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Projects/biopax-pattern/All-Data.owl"));
//		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Projects/biopax-pattern/Captured-state-change.owl"));

		PathAccessor paPart = new PathAccessor("SmallMolecule/participantOf:Conversion");
		PathAccessor paPartCA = new PathAccessor("SmallMolecule/participantOf:ComplexAssembly");
		PathAccessor paEff = new PathAccessor ("SmallMolecule/controllerOf/controlled*:Conversion");
		PathAccessor paComp = new PathAccessor ("SmallMolecule/componentOf:Complex");

		List<SMHolder> list = new ArrayList<SMHolder>();
		for (SmallMolecule smr : model.getObjects(SmallMolecule.class))
		{
			list.add(new SMHolder(smr,
				paPart.getValueFromBean(smr).size() - paPartCA.getValueFromBean(smr).size(),
				paEff.getValueFromBean(smr).size(),
				paComp.getValueFromBean(smr).size()));
		}
		Collections.sort(list);

		for (SMHolder s : list)
		{
			if (s.partDegree > 5) System.out.println(s);
		}
	}

	class SMHolder implements Comparable
	{
		Named smr;
		int partDegree;
		int effDegree;
		int compDegree;

		SMHolder(Named smr, int partDegree, int effDegree, int compDegree)
		{
			this.smr = smr;
			this.partDegree = partDegree;
			this.effDegree = effDegree;
			this.compDegree = compDegree;
		}

		@Override
		public int compareTo(Object o)
		{
			if (o instanceof SMHolder)
			{
				SMHolder s = (SMHolder) o;
				return - (partDegree) + (s.partDegree);
			}
			return 0;
		}

		@Override
		public String toString()
		{
			return partDegree + "  \t" + effDegree + "  \t" + compDegree + "\t" + smr.getDisplayName();
		}

		int degree()
		{
			return partDegree + effDegree;
		}
	}

	@Test
	@Ignore
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
