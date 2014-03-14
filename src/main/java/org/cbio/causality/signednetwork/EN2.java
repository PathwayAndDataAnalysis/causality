package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.miner.ControlsExpressionWithConvMiner;

/**
 * @author Ozgun Babur
 */
public class EN2 extends ControlsExpressionWithConvMiner
{
	public EN2()
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
