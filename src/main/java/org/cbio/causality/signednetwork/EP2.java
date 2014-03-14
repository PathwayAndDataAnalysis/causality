package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.miner.ControlsExpressionWithConvMiner;

/**
 * @author Ozgun Babur
 */
public class EP2 extends ControlsExpressionWithConvMiner
{
	public EP2()
	{
		setType(SignedType.UPREGULATES_EXPRESSION);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), "Control");
		return p;
	}
}
