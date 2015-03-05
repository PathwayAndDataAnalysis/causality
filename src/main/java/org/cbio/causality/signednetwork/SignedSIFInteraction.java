package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SignedSIFInteraction extends SIFInteraction
{
	Set<String> changedPhospho;
	Map<String, Set<BioPAXElement>> p2med;

	public SignedSIFInteraction(String sourceID, String targetID, BioPAXElement sourceER, BioPAXElement targetER, SIFType type,
		Set<BioPAXElement> mediators, Set<BioPAXElement> sourcePEs, Set<BioPAXElement> targetPEs,
		Set<String> changedPhospho)
	{
		super(sourceID, targetID, sourceER, targetER, type, mediators, sourcePEs, targetPEs);
		this.changedPhospho = changedPhospho;
		this.p2med = new HashMap<String, Set<BioPAXElement>>();
		for (String ph : changedPhospho)
		{
			p2med.put(ph, new HashSet<BioPAXElement>(mediators));
		}
	}

	@Override
	public void mergeWith(SIFInteraction equivalent)
	{
		if (equivalent instanceof SignedSIFInteraction)
		{
			SignedSIFInteraction ss = (SignedSIFInteraction) equivalent;
			super.mergeWith(equivalent);
			this.changedPhospho.addAll(ss.changedPhospho);
			for (String ph : ss.p2med.keySet())
			{
				if (p2med.containsKey(ph)) p2med.get(ph).addAll(ss.p2med.get(ph));
				else p2med.put(ph, ss.p2med.get(ph));
			}
		}
	}

	@Override
	public String toString(boolean withMediators)
	{
		String s = super.toString(withMediators) + "\t";

		for (String site : changedPhospho)
		{
			s += site + ";";
		}
		if (s.endsWith(";")) s = s.substring(0, s.length() - 1) + "\t";

		for (String pmid : getPubmedIDs())
		{
			s += pmid + " ";
		}
		if (s.endsWith(" ")) s = s.substring(0, s.length() - 1);

		return s;
	}
}
