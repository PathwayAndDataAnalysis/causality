package org.cbio.causality.network;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public interface InteractionProvider
{
	public Set<String> getInteractions(String symbol);
}
