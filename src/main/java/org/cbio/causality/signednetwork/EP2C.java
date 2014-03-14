package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.NOT;
import org.biopax.paxtools.pattern.constraint.Type;
import org.biopax.paxtools.pattern.miner.ControlsExpressionWithConvMiner;

/**
 * @author Ozgun Babur
 */
public class EP2C extends EP2
{
	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(new Type(Complex.class)), "TF PE");
		return p;
	}
}
