package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SignedSIFInteraction extends SIFInteraction
{
	Set<String> changedPhospho;
	public SignedSIFInteraction(BioPAXElement sourceER, BioPAXElement targetER, SIFType type,
		Set<BioPAXElement> mediators, Set<BioPAXElement> sourcePEs, Set<BioPAXElement> targetPEs,
		IDFetcher fetcher, Set<String> changedPhospho)
	{
		super(sourceER, targetER, type, mediators, sourcePEs, targetPEs, fetcher);
		this.changedPhospho = changedPhospho;
	}

	@Override
	public void mergeWith(SIFInteraction equivalent)
	{
		super.mergeWith(equivalent);
		this.changedPhospho.addAll(((SignedSIFInteraction) equivalent).changedPhospho);
	}

	@Override
	public String toString(boolean withMediators)
	{
		String s = super.toString(withMediators) + "\t";

		for (String site : changedPhospho)
		{
			s += site + ";";
		}
		s = s.substring(0, s.length() - 1);
		return s;
	}
}
