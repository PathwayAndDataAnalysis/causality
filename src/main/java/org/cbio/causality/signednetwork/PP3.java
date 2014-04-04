package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.ConBox;
import org.biopax.paxtools.pattern.constraint.ModificationChangeConstraint;
import org.biopax.paxtools.pattern.constraint.NOT;
import org.biopax.paxtools.pattern.miner.CSCOBothControllerAndParticipantMiner;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.HashSet;

/**
 * @author Ozgun Babur
 */
public class PP3 extends CSCOBothControllerAndParticipantMiner
{
	public PP3()
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

	@Override
	public SIFInteraction createSIFInteraction(Match m, IDFetcher fetcher)
	{
		return new SignedSIFInteraction(m.get(this.getSourceLabel(), getPattern()),
			m.get(this.getTargetLabel(), getPattern()), this.getSIFType(),
			new HashSet<BioPAXElement>(m.get(getMediatorLabels(), getPattern())),
			new HashSet<BioPAXElement>(m.get(getSourcePELabels(), getPattern())),
			new HashSet<BioPAXElement>(m.get(getTargetPELabels(), getPattern())),
			fetcher, DifferentialModificationUtil.collectChangedPhosphorylationSites(
				(PhysicalEntity) m.get("input simple PE", getPattern()),
				(PhysicalEntity) m.get("output simple PE", getPattern()),
				getSIFType().equals(SignedType.PHOSPHORYLATES)));
	}
}
