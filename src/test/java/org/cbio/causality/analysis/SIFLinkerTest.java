package org.cbio.causality.analysis;

import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SIFLinkerTest
{
	@Test
	public void testLinking()
	{
		SIFLinker linker = new SIFLinker();
		linker.load("/home/ozgun/Projects/mutex/src/main/resources/network.txt");

		Set<String> from = new HashSet<String>(Arrays.asList("TP53", "PTEN", "CTNNB1", "PIK3CA", "ARID1A", "KRAS", "ARID5A"));

		List<String> set = linker.link(from, from, 0);

		for (String s : set)
		{
			System.out.println(s);
		}
	}
}
