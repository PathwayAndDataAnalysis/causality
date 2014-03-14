package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.ConBox;
import org.biopax.paxtools.pattern.constraint.ModificationChangeConstraint;
import org.biopax.paxtools.pattern.constraint.NOT;
import org.biopax.paxtools.pattern.miner.CSCOButIsParticipantMiner;

/**
 * @author Ozgun Babur
 */
public class PP2 extends CSCOButIsParticipantMiner
{
	public PP2()
	{
		setType(SignedType.PHOSPHORYLATES);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(ConBox.linkToSimple()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSimple()), "output PE", "input simple PE");
		p.add(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"),
			"input simple PE", "output simple PE");
		return p;
	}
}
