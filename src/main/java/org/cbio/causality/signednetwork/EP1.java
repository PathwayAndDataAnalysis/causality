package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.AND;
import org.biopax.paxtools.pattern.constraint.ModificationChangeConstraint;
import org.biopax.paxtools.pattern.constraint.OR;
import org.biopax.paxtools.pattern.miner.ControlsExpressionMiner;
import org.biopax.paxtools.pattern.miner.ControlsStateChangeOfMiner;

/**
 * @author Ozgun Babur
 */
public class EP1 extends ControlsExpressionMiner
{
	public EP1()
	{
		setType(SignedType.UPREGULATES_EXPRESSION);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NotLessDirect(), "TF ER", "TF SPE", "TF PE", "TempReac");
		p.add(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), "Control");
		return p;
	}
}
