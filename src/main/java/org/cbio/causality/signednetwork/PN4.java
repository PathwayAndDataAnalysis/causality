package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.CSCOThroughControllingSmallMoleculeMiner;

/**
 * @author Ozgun Babur
 */
public class PN4 extends PP4
{
	public PN4()
	{
		setType(SignedType.DEPHOSPHORYLATES);
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
					new MappedConst(
						new OR(
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 1)
								), 0, 1),
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 1)
								), 0, 1)
						), 0, 1
					),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"), 2, 3)
				), 0, 1, 2, 3),
			new MappedConst(
				new AND(
					new MappedConst(
						new OR(
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 1)
								), 0, 1),
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 1)
								), 0, 1)
						), 0, 1
					),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"), 2, 3)
				), 0, 1, 2, 3)
			), "upper Control", "Control", "input simple PE", "output simple PE");
		return p;
	}
}
