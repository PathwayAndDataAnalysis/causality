package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class GraphTest
{
	@Test
	@Ignore
	public void printDegreeDistribution() throws FileNotFoundException
	{
		String file = "/home/ozgun/Desktop/PC.sif";
		Graph trav = new Graph();

		for (SIFType type : SIFEnum.values())
		{
			System.out.println("\ntype = " + type.getTag());
			trav.load(file,
				type.isDirected() ? Collections.<String>emptySet() : Collections.singleton(type.getTag()),
				type.isDirected() ? Collections.singleton(type.getTag()) : Collections.<String>emptySet());

			System.out.println("\nindegree");
			printDegrees(trav.getDegreeDistibution(true));
			System.out.println("\noutdegree");
			printDegrees(trav.getDegreeDistibution(false));
			System.out.println("\ntotal");
			printDegrees(trav.getDegreeDistibution());

			trav.clear();
		}
	}

	private void printDegrees(Map<Integer, Integer> cnts)
	{
		List<Integer> degrees = new ArrayList<Integer>(cnts.keySet());
		Collections.sort(degrees);
		int totalDegrees = 0;
		int totalCnt = 0;
		for (Integer degree : degrees)
		{
			Integer cnt = cnts.get(degree);
			System.out.println(degree + "\t" + cnt);
			totalDegrees += degree * cnt;
			totalCnt += cnt;
		}

		System.out.println("Average Degree = " + (totalDegrees / (double) totalCnt));
	}
}
