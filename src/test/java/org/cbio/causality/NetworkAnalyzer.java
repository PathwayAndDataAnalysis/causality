package org.cbio.causality;

import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.Kronometre;
import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.*;
import org.biopax.paxtools.pattern.constraint.*;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * This is a temporary class for analyzing the big picture we have. I should have made this an 
 * independent project and use paxtools, but I am lazy. 
 * 
 * @author Ozgun Babur
 */
public class NetworkAnalyzer
{
	@Test
	@Ignore
	public void lookForGeneralInhibitors() throws FileNotFoundException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Desktop/all.owl"));

		for (ProteinReference er : model.getObjects(ProteinReference.class))
		{
			PathAccessor pa = new PathAccessor("EntityReference/entityReferenceOf/componentOf*");
			Set comps = pa.getValueFromBean(er);

			List<Match> matches = Searcher.search(er, PatternBox.inSameComplex());

			Set<EntityReference> interacting = collect(matches, 4);
			
			matches = Searcher.search(er, PatternBox.inSameComplexHavingTransActivity());

			Set<EntityReference> transAct = collect(matches, 4);
			
			matches = Searcher.search(er, PatternBox.inSameComplexEffectingConversion());

			Set<EntityReference> effConv = collect(matches, 4);

			//todo revise this test. below method is gone from patternbox
//			matches = Searcher.search(er, PatternBox.hasActivity());

			Set<Control> controls = collect(matches, 2);

			System.out.println(comps.size() + "\t" + interacting.size() + "\t" + transAct.size() + "\t" + effConv.size() + "\t" + controls.size() + "\t" + er.getDisplayName());
		}
	}
	
	protected Set collect(Collection<Match> matches, int index)
	{
		Set set = new HashSet();
		for (Match match : matches)
		{
			set.add(match.get(index));
		}
		return set;
	}

	@Test
	@Ignore
	public void clip() throws FileNotFoundException
	{
		Kronometre k = new Kronometre();
		
		Pattern p = new Pattern(null, ""); // pattern here

		Searcher.searchInFile(p, "/home/ozgun/Desktop/PC.owl", "/home/ozgun/Desktop/pattern-matches/ChangedToEffect.owl", 100, 2);

		k.stop();
		k.print();
		
//		SimpleIOHandler h = new SimpleIOHandler();
//		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Desktop/cpath2_prepared.owl"));
//
//		List<Match> mathes = Searcher.search(model.getByID("urn:miriam:uniprot:Q01094"), p);
//		System.out.println("mathes.size() = " + mathes.size());
	}


	@Test
	@Ignore
	public void createSIF() throws IOException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Desktop/cpath2.owl"));

		Pattern bindPattern = PatternBox.bindsTo();

		Pattern stChPattern = PatternBox.controlsStateChange();
		stChPattern.add(new Type(Protein.class), "controller ER");
		stChPattern.add(new Type(Protein.class), "changed ER");

		Pattern ppiPattern = PatternBox.molecularInteraction();
		Pattern trConvPattern = PatternBox.controlsExpressionWithConversion();
		Pattern trTempPattern = PatternBox.controlsExpressionWithTemplateReac();
		Pattern degredPattern = PatternBox.controlsStateChangeThroughDegradation();

		BufferedWriter writer = new BufferedWriter(new FileWriter("/home/ozgun/Desktop/SIF.txt"));

		int i = 0;
		for (ProteinReference pr : model.getObjects(ProteinReference.class))
		{
			if (++i % 10 == 0) System.out.print(".");
			if (i % 1000 == 0) System.out.println();

			if (getSymbol(pr) == null) continue;
			
			List<Match> matches = Searcher.search(pr, bindPattern);
			matches.addAll(Searcher.search(pr, ppiPattern));
			recordSIF(matches, "BINDS_TO", writer);

			matches = Searcher.search(pr, stChPattern);
			recordSIF(matches, "STATE_CHANGE", writer);

			matches = Searcher.search(pr, trConvPattern);
			matches.addAll(Searcher.search(pr, trTempPattern));
			recordSIF(matches, "TRANSCRIPTION", writer);

			matches = Searcher.search(pr, degredPattern);
			recordSIF(matches, "DEGRADATION", writer);
		}
		
		writer.close();
	}

	private void recordSIF(List<Match> matchMap, String type, BufferedWriter writer) 
		throws IOException
	{
		Set<String> existing = new HashSet<String>();

		for (Match match : matchMap)
		{
			ProteinReference pr = (ProteinReference) match.getFirst();
			String s1 = getSymbol(pr);

			EntityReference er = (EntityReference) match.getLast();
			if (er instanceof ProteinReference)
			{
				String s2 = getSymbol((ProteinReference) er);
				if (s2 != null)
				{
					String line = s1 + "\t" + type + "\t" + s2;
					
					if (!existing.contains(line))
					{
						writer.write(line + "\n");
						existing.add(line);
					}
				}
			}
		}
	}

	private String getSymbol(ProteinReference pr)
	{
		for (Xref xref : pr.getXref())
		{
			if (xref.getDb().equals("HGNC"))
			{
				String id = xref.getId().substring(xref.getId().indexOf(":") + 1);
				return HGNC.getSymbol(id);
			}
		}
		return null;
	}

}
