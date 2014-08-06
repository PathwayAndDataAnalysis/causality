package org.cbio.causality.util;

/**
 * A base directory provider class for the resources to be downloaded.
 *
 * @author Ozgun Babur
 */
public class BaseDir
{
	private static String dir = "";

	public static String getDir()
	{
		return dir;
	}

	public static void setDir(String dir)
	{
		BaseDir.dir = dir;
	}
}
