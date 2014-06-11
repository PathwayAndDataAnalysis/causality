package org.cbio.causality.analysis;

import org.cbio.causality.model.RPPAData;
import org.cbio.causality.network.SignedPC;
import org.cbio.causality.signednetwork.SignedType;
import org.cbio.causality.util.CollectionUtil;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RPPANetworkMapper
{
	public static List<Relation> map(Collection<RPPAData> data)
	{
		List<Relation> rels = new ArrayList<Relation>();
		Map<SignedType,Graph> graph = SignedPC.getAllGraphs();

		for (RPPAData d1 : data)
		{
			for (RPPAData d2 : data)
			{
				if (d1 == d2 || d1.genes.equals(d2.genes)) continue;

				for (String source : d1.genes)
				{
					for (String target : d2.genes)
					{
						if (source.equals(target)) continue;

						for (SignedType type : graph.keySet())
						{
							if (graph.get(type).getDownstream(source).contains(target))
							{
								boolean match = false;
								int sign = 0;

								switch (type)
								{
									case PHOSPHORYLATES:
									{
										if (d2.isPhospho())
										{
											match = true;
											sign = d1.getSelfEffect();
										}
									} break;
									case DEPHOSPHORYLATES:
									{
										if (d2.isPhospho())
										{
											match = true;
											sign = -d1.getSelfEffect();
										}
									} break;
									case UPREGULATES_EXPRESSION:
									{
										if (!d2.isPhospho())
										{
											match = true;
											sign = d1.getSelfEffect();
										}
									} break;
									case DOWNREGULATES_EXPRESSION:
									{
										if (!d2.isPhospho())
										{
											match = true;
											sign = -d1.getSelfEffect();
										}
									} break;
								}

								if (match)
								{
									Set<String> sites = (graph.get(type) instanceof PhosphoGraph) ?
										((PhosphoGraph) graph.get(type)).getSites(source, target) :
										null;

									rels.add(new Relation(source, target, type, d1, d2, sign,
										graph.get(type).getMediatorsInString(source, target),
										sites));
								}
							}
						}
					}
				}
			}
		}
		return rels;
	}

	public static void removeConflicting(List<Relation> rels)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (rel.dataChangesAsUnxpected()) iter.remove();
		}
	}

	public static void keepMatching(List<Relation> rels)
	{
		Iterator<Relation> iter = rels.iterator();
		while (iter.hasNext())
		{
			Relation rel = iter.next();
			if (!rel.dataChangesAsExpected()) iter.remove();
		}
	}
}
