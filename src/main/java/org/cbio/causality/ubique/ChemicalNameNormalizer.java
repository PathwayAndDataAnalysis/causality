package org.cbio.causality.ubique;

/**
 * @author Ozgun Babur
 */
public class ChemicalNameNormalizer
{
	public static String nornmalize(String name)
	{
		if (name == null) return null;

		name = name.toLowerCase();
		if (name.endsWith(")") && name.lastIndexOf("(") > name.length() - 6)
			name = name.substring(0, name.lastIndexOf("(")).trim();
		if (name.startsWith("l-"))
			name = name.substring(2).trim();
		if (name.endsWith("ic acid"))
			name = name.substring(0, name.length()-7) + "ate";

		return name;
	}
}
