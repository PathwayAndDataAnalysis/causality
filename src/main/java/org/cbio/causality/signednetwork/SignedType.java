package org.cbio.causality.signednetwork;

import org.biopax.paxtools.pattern.miner.SIFMiner;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.util.Arrays;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public enum SignedType implements SIFType
{
	PHOSPHORYLATES("First protein positively affects phosphorylation of the second protein.",
		true, true),
	DEPHOSPHORYLATES("First protein negatively affects phosphorylation of the second protein.",
		true, true),
	UPREGULATES_EXPRESSION("First protein positively affects expression of the second protein.",
		true, false),
	DOWNREGULATES_EXPRESSION("First protein negatively affects expression of the second protein.",
		true, false),
	;

	/**
	 * Constructor with parameters.
	 * @param description description of the edge type
	 * @param directed whether the edge type is directed
	 */
	private SignedType(String description, boolean directed, boolean phospho,
		Class<? extends SIFMiner>... miners)
	{
		this.description = description;
		this.directed = directed;
		this.phospho = phospho;
		this.miners = Arrays.asList(miners);
	}

	/**
	 * Description of the SIF type.
	 */
	private String description;

	/**
	 * Some SIF edges are directed and others are not.
	 */
	private boolean directed;

	/**
	 * Whether this relation is related to a phosphorylation event.
	 */
	private boolean phospho;

	/**
	 * SIF Miners to use during a search.
	 */
	private List<Class<? extends SIFMiner>> miners;

	/**
	 * Tag of a SIF type is derived from the enum name.
	 * @return tag
	 */
	public String getTag()
	{
		return name().toLowerCase().replaceAll("_", "-");
	}

	/**
	 * Asks if the edge is directed.
	 * @return true if directed
	 */
	public boolean isDirected()
	{
		return directed;
	}

	public boolean isPhospho()
	{
		return phospho;
	}

	/**
	 * Gets the description of the SIF type.
	 * @return description
	 */
	public String getDescription()
	{
		return description;
	}

	@Override
	public List<Class<? extends SIFMiner>> getMiners()
	{
		return miners;
	}

	public static SignedType typeOf(String tag)
	{
		tag = tag.toUpperCase().replaceAll("-", "_");
		SignedType type = null;
		try
		{
			type = valueOf(tag);
		}
		catch (IllegalArgumentException e){}
		return type;
	}
}
