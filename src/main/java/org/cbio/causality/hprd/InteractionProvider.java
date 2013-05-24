package org.cbio.causality.hprd;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public interface InteractionProvider
{
	public Set<String> getInteractions(String symbol);
}
