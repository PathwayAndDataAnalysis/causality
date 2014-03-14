package org.cbio.causality.signednetwork;

import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.constraint.ConstraintAdapter;

/**
 * @author Ozgun Babur
 */
public class ControlSignConstraint extends ConstraintAdapter
{
	protected Sign desiredSign;

	public ControlSignConstraint(Sign desiredSign)
	{
		super(1);
		this.desiredSign = desiredSign;
	}

	@Override
	public boolean satisfies(Match match, int... ind)
	{
		Control ctrl = (Control) match.get(ind[0]);
		boolean neg = ctrl.getControlType() != null &&
			ctrl.getControlType().toString().startsWith("I");

		switch (desiredSign)
		{
			case POSITIVE: return !neg;
			case NEGATIVE: return neg;
		}
		return false;
	}

	public enum Sign
	{
		POSITIVE,
		NEGATIVE
	}
}
