package org.cbio.causality.analysis;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Named;
import org.biopax.paxtools.model.level3.SmallMolecule;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.cbio.causality.util.Histogram;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class TraverseTest
{
	@Test
	@Ignore
	public void printDegreeDistribution() throws FileNotFoundException
	{
		Traverse trav = new Traverse();
		trav.load("SIF.txt", Collections.EMPTY_SET, new HashSet<String>(Arrays.asList(
			SIFType.CONTROLS_STATE_CHANGE.getTag(),
			SIFType.CONTROLS_EXPRESSION.getTag(),
			SIFType.CONTROLS_DEGRADATION.getTag())));

		Histogram hin = new Histogram(5);
		Histogram hout = new Histogram(5);

		for (String s : trav.getSymbols())
		{
			Set<String> upstr = trav.getUpstream(s);
			int indegree = upstr.size();
			Set<String> dwstr = trav.getDownstream(s);
			int outdegree = dwstr.size();
			upstr.retainAll(dwstr);
			int comm = upstr.size();

			hin.count(indegree);
			hout.count(outdegree);

			if (indegree + outdegree > 100)
			{
				System.out.println(s + "\t" + indegree + "\t" + outdegree + "\t" + comm);
			}
		}

//		hin.printTogether(hout);
	}
}
