package org.cbio.causality.signednetwork;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFSearcher;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.cbio.causality.util.Kronometre;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class Generator
{
	public static void main(String[] args) throws Throwable
	{
		Kronometre kron = new Kronometre();

		SimpleIOHandler h = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = h.convertFromOWL(new FileInputStream("../biopax-pattern/All-Data.owl"));
//		Model model = h.convertFromOWL(new FileInputStream("../biopax-pattern/PANTHER.owl"));
		Blacklist blacklist = new Blacklist("../biopax-pattern/blacklist.txt");
		SIFSearcher searcher;

		// prepare phospho-graph

		searcher = new SIFSearcher(new PP1(), new PP2(), new PP3(), new PP4());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> pp = searcher.searchSIF(model);

		searcher = new SIFSearcher(new PP1C(), new PP4C());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> ppc = searcher.searchSIF(model);

		searcher = new SIFSearcher(new PN1(), new PN2(), new PN3(), new PN4());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> pn = searcher.searchSIF(model);

		searcher = new SIFSearcher(new PN1C(), new PN4C());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> pnc = searcher.searchSIF(model);

		decideConflicting(pp, pn, ppc, pnc);

		// prepare expression graph

		searcher = new SIFSearcher(new EP1(), new EP2());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> ep = searcher.searchSIF(model);

		searcher = new SIFSearcher(new EP1C(), new EP2C());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> epc = searcher.searchSIF(model);

		searcher = new SIFSearcher(new EN1(), new EN2());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> en = searcher.searchSIF(model);

		searcher = new SIFSearcher(new EN1C(), new EN2C());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> enc = searcher.searchSIF(model);

		decideConflicting(ep, en, epc, enc);

		Set<SIFInteraction> sifs = new HashSet<SIFInteraction>();
		sifs.addAll(pp);
		sifs.addAll(pn);
		sifs.addAll(ep);
		sifs.addAll(en);

		BufferedWriter writer = new BufferedWriter(new FileWriter("SignedPC.sif"));

		for (SIFInteraction sif : sifs)
		{
			writer.write(sif.toString(true) + "\n");
		}

		writer.close();
		kron.stop();
		kron.print();
	}

	private static void decideConflicting(Set<SIFInteraction> pos, Set<SIFInteraction> neg,
		Set<SIFInteraction> posC, Set<SIFInteraction> negC)
	{
		Set<String> ov = getWOTypes(pos);
		ov.retainAll(getWOTypes(neg));
		System.out.println("\noverlap size = " + ov.size());

		Set<String> pc = getWOTypes(posC);
		Set<String> nc = getWOTypes(negC);

		Set<String> tmp = new HashSet<String>(pc);
		pc.removeAll(nc);
		nc.removeAll(tmp);

		Map<String, Boolean> decision = new HashMap<String, Boolean>();

		for (String s : ov)
		{
			if (pc.contains(s)) decision.put(s, true);
			else if (nc.contains(s)) decision.put(s, false);
		}

		System.out.println("decision.size() = " + decision.size());
		for (String s : decision.keySet())
		{
			System.out.println(s + "\t" + decision.get(s));
		}

		for (SIFInteraction sif : new HashSet<SIFInteraction>(pos))
		{
			String key = woType(sif);
			if (decision.containsKey(key) && !decision.get(key))
			{
				pos.remove(sif);
			}
		}
		for (SIFInteraction sif : new HashSet<SIFInteraction>(neg))
		{
			String key = woType(sif);
			if (decision.containsKey(key) && decision.get(key))
			{
				neg.remove(sif);
			}
		}
	}

	private static String woType(SIFInteraction sif)
	{
		return sif.sourceID + "\t" + sif.targetID;
	}

	private static Set<String> getWOTypes(Set<SIFInteraction> sifs)
	{
		Set<String> set = new HashSet<String>(sifs.size());

		for (SIFInteraction sif : sifs)
		{
			set.add(woType(sif));
		}
		return set;
	}
}
