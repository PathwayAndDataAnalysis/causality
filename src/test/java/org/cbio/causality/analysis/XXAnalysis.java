package org.cbio.causality.analysis;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.pattern.miner.*;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.biopax.paxtools.trove.TProvider;
import org.biopax.paxtools.util.BPCollections;
import org.cbio.causality.util.Kronometre;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

/**
 * @author Ozgun Babur
 */
public class XXAnalysis
{
	@Test
	@Ignore
	public void doit() throws FileNotFoundException
	{
		Kronometre kron = new Kronometre();

		BPCollections.I.setProvider(new TProvider());

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream
			("/home/ozgun/Projects/biopax-pattern/All-Data.owl"));
//			("/home/ozgun/Projects/biopax-pattern/Reactome.owl"));

		Blacklist blacklist = new Blacklist("/home/ozgun/Projects/biopax-pattern/blacklist.txt");

		SIFSearcher s = new SIFSearcher(SIFEnum.values());
//		SIFSearcher s = new SIFSearcher(SIFEnum.CONTROLS_TRANSPORT_OF);
		s.setBlacklist(blacklist);

		s.searchSIF(model, new FileOutputStream("SIFWithLoc.sif"),
			new CustomFormat(
				OutputColumn.Type.SOURCE_LOC.name(),
				OutputColumn.Type.TARGET_LOC.name(),
				OutputColumn.Type.RESOURCE.name(),
				OutputColumn.Type.PATHWAY.name()
				));

		kron.stop();
		kron.print();
		Assert.assertTrue(2 + 2 == 4);
	}
}
