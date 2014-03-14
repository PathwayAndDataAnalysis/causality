package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.miner.ControlsExpressionMiner;

/**
 * @author Ozgun Babur
 */
public class EN1 extends ControlsExpressionMiner
{
	public EN1()
	{
		setType(SignedType.DOWNREGULATES_EXPRESSION);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), "Control");
		return p;
	}
}
