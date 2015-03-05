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
import org.biopax.paxtools.pattern.miner.SIFType;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.HashSet;
import java.util.Set;

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
		p.add(new NOT(ConBox.linkToSpecific()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSpecific()), "output PE", "input simple PE");
		p.add(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"),
			"input simple PE", "output simple PE");
		return p;
	}

	@Override
	public Set<SIFInteraction> createSIFInteraction(Match m, IDFetcher fetcher)
	{
		BioPAXElement sourceER = m.get(this.getSourceLabel(), getPattern());
		BioPAXElement targetER = m.get(this.getTargetLabel(), getPattern());

		Set<String> sources = fetcher.fetchID(sourceER);
		Set<String> targets = fetcher.fetchID(targetER);

		SIFType sifType = this.getSIFType();

		Set<SIFInteraction> set = new HashSet<SIFInteraction>();

		for (String source : sources)
		{
			for (String target : targets)
			{
				if (source.equals(target)) continue;

				set.add(new SignedSIFInteraction(source, target, sourceER, targetER, sifType,
					new HashSet<BioPAXElement>(m.get(getMediatorLabels(), getPattern())),
					new HashSet<BioPAXElement>(m.get(getSourcePELabels(), getPattern())),
					new HashSet<BioPAXElement>(m.get(getTargetPELabels(), getPattern())),
					DifferentialModificationUtil.collectChangedPhosphorylationSites(
						(PhysicalEntity) m.get("input simple PE", getPattern()),
						(PhysicalEntity) m.get("output simple PE", getPattern()),
						getSIFType().equals(SignedType.PHOSPHORYLATES))));
			}
		}
		return set;
	}
}
