package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.network.MSigDBTFT;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.network.SignaLink;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
@Ignore
public class EdsAnalysis
{
	@Test
	@Ignore
	public void findGraphAroundSmallMol() throws IOException
	{
//		Set<String> smNames = new HashSet<String>(Arrays.asList("sorbitol", "glucitol"));
		Set<String> smNames = new HashSet<String>(Arrays.asList("d-fructose"));

		Set<String> relations = new HashSet<String>();

		Graph ccbGraph = PathwayCommons.getGraph(SIFEnum.CONSUMPTION_CONTROLLED_BY);
		Graph cpoGraph = PathwayCommons.getGraph(SIFEnum.CONTROLS_PRODUCTION_OF);

		cpoGraph.putRelation("AKR1A1", "D-glucitol", true);

		Set<String> sms = findNodes(ccbGraph, smNames);
		sms.addAll(findNodes(cpoGraph, smNames));

		System.out.println("sms = " + sms);

		int d0 = 3;

		Set<String> prod = traverseSMP(sms, ccbGraph, cpoGraph, false, relations, d0);
		System.out.println("prod.size() = " + prod.size());

		Set<String> cons = traverseSMP(sms, ccbGraph, cpoGraph, true, relations, d0);
		System.out.println("cons.size() = " + cons.size());

		Graph catPreGraph = PathwayCommons.getGraph(SIFEnum.CATALYSIS_PRECEDES);

		int d1 = 0;

		Set<String> prodFurther = catPreGraph.getUpstream(prod, d1);
		System.out.println("prodFurther.size() = " + prodFurther.size());
		Set<String> consFurther = catPreGraph.getDownstream(cons, d1);
		System.out.println("consFurther.size() = " + consFurther.size());

		Set<String> met = new HashSet<String>(prod);
		met.addAll(cons);
		met.addAll(prodFurther);
		met.addAll(consFurther);
		System.out.println("met.size() = " + met.size());

		if (d1 > 0) relations.addAll(catPreGraph.toString(met, met));

		Graph stchGraph = PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
		stchGraph.merge(SPIKE.getGraphPostTl());
		stchGraph.merge(SignaLink.getGraphPostTl());

		int d2 = 0;

		Set<String> regs = stchGraph.getUpstream(met, d2);
		System.out.println("regs.size() = " + regs.size());

		if (d2 > 1) relations.addAll(stchGraph.toString(regs, regs));
		relations.addAll(stchGraph.toString(regs, met));

		Set<String> metAndRegs = new HashSet<String>(met);
		metAndRegs.addAll(regs);
		System.out.println("metAndRegs.size() = " + metAndRegs.size());

		Graph expGraph = PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		expGraph.merge(SPIKE.getGraphTR());
		expGraph.merge(SignaLink.getGraphTR());
		expGraph.merge(MSigDBTFT.getGraph());

		int d3 = 0;

		Set<String> tfs = expGraph.getUpstream(metAndRegs, d3);
		System.out.println("tfs.size() = " + tfs.size());

		if (d3 > 1) relations.addAll(expGraph.toString(tfs, tfs));
		relations.addAll(expGraph.toString(tfs, metAndRegs));

		BufferedWriter writer = new BufferedWriter(new FileWriter("eds-network.sif"));
		for (String relation : relations)
		{
			writer.write(relation + "\n");
		}
		writer.close();
	}

	private Set<String> traverseSMP(Set<String> query, Graph ccbGraph, Graph cpoGraph, boolean upstream,
		Set<String> relations, int depth)
	{
		Set<String> neighborProts = new HashSet<String>();
		Set<String> neighborChems = new HashSet<String>();
		Set<String> currentQ = new HashSet<String>(query);
		for (int i = 0; i < depth; i++)
		{
			if (upstream)
			{
				Set<String> proteins = cpoGraph.getUpstream(currentQ);
				Set<String> chems = ccbGraph.getUpstream(proteins);
				proteins.removeAll(neighborProts);
				chems.removeAll(neighborChems);

				relations.addAll(cpoGraph.toString(proteins, currentQ));
				relations.addAll(ccbGraph.toString(chems, proteins));

				neighborProts.addAll(proteins);
				neighborChems.addAll(chems);
				currentQ.clear();
				currentQ.addAll(chems);
			}
			else
			{
				Set<String> proteins = ccbGraph.getDownstream(currentQ);
				Set<String> chems = cpoGraph.getDownstream(proteins);
				proteins.removeAll(neighborProts);
				chems.removeAll(neighborChems);

				relations.addAll(ccbGraph.toString(currentQ, proteins));
				relations.addAll(cpoGraph.toString(proteins, chems));

				neighborProts.addAll(proteins);
				neighborChems.addAll(chems);
				currentQ.clear();
				currentQ.addAll(chems);
			}
		}
		return neighborProts;
	}

	private Set<String> findNodes(Graph graph, Set<String> query)
	{
		Set<String> result = new HashSet<String>();

		for (String s : graph.getSymbols())
		{
			for (String q : query)
			{
//				if (s.toLowerCase().contains(q))
				if (s.toLowerCase().equals(q))
				{
					result.add(s);
					break;
				}
			}
		}
		return result;
	}
}
