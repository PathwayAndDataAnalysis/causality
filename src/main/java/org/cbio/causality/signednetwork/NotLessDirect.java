package org.cbio.causality.signednetwork;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.model.level3.EntityReference;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.constraint.ConstraintAdapter;

/**
 * Makes sure that the ER is not controlling the given interaction more directly through another
 * state.
 *
 * var0: Upstream ER
 * var1: Upstream simple PE
 * var2: Upstream complex PE
 * var3: The controlled conversion or template reaction
 *
 * @author Ozgun Babur
 */
public class NotLessDirect extends ConstraintAdapter
{
	private final static PathAccessor pa = new PathAccessor(
		"EntityReference/entityReferenceOf/controllerOf/controlled");

	public NotLessDirect()
	{
		super(4);
	}

	@Override
	public boolean satisfies(Match match, int... ind)
	{
		EntityReference er = (EntityReference) match.get(ind[0]);
		PhysicalEntity spe = (PhysicalEntity) match.get(ind[1]);
		PhysicalEntity pe = (PhysicalEntity) match.get(ind[2]);
		Interaction inter = (Interaction) match.get(ind[3]);

		return pe == spe || !pa.getValueFromBean(er).contains(inter);
	}
}
