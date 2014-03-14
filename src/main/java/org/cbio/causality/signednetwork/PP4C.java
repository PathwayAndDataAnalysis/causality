package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.CSCOThroughControllingSmallMoleculeMiner;

/**
 * @author Ozgun Babur
 */
public class PP4C extends PP4
{
	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(new Type(Complex.class)), "upper controller PE");
		return p;
	}
}
