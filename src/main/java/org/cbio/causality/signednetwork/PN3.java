package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.ConBox;
import org.biopax.paxtools.pattern.constraint.ModificationChangeConstraint;
import org.biopax.paxtools.pattern.constraint.NOT;
import org.biopax.paxtools.pattern.miner.CSCOBothControllerAndParticipantMiner;

/**
 * @author Ozgun Babur
 */
public class PN3 extends PP3
{
	public PN3()
	{
		setType(SignedType.DEPHOSPHORYLATES);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(ConBox.linkToSimple()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSimple()), "output PE", "input simple PE");
		p.add(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"),
			"input simple PE", "output simple PE");
		return p;
	}
}
