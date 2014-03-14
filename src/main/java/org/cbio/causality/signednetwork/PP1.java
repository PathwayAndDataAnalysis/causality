package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.ControlsStateChangeOfMiner;

/**
 * @author Ozgun Babur
 */
public class PP1 extends ControlsStateChangeOfMiner
{
	public PP1()
	{
		setType(SignedType.PHOSPHORYLATES);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(ConBox.linkToSimple()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSimple()), "output PE", "input simple PE");
		p.add(new OR(
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"), 1, 2)
				), 0, 1, 2),
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"), 1, 2)
				), 0, 1, 2)
			), "Control", "input simple PE", "output simple PE");
		return p;
	}
}
