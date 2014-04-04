package org.cbio.causality.analysis;

import org.cbio.causality.model.RPPAData;
import org.cbio.causality.network.SignedPC;
import org.cbio.causality.signednetwork.SignedType;

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

				for (String g1 : d1.genes)
				{
					for (String g2 : d2.genes)
					{
						assert !g1.equals(g2);

						for (SignedType type : graph.keySet())
						{
							if (graph.get(type).getDownstream(g1).contains(g2))
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
									rels.add(new Relation(g1, g2, type, d1, d2, sign,
										graph.get(type).getMediatorsInString(g1, g2)));
								}
							}
						}
					}
				}
			}
		}
		return rels;
	}


	public static class Relation
	{
		public Relation(String source, String target, SignedType edgeType,
			RPPAData sourceData, RPPAData targetData, int corrSign, String mediators)
		{
			this.source = source;
			this.target = target;
			this.edgeType = edgeType;
			this.sourceData = sourceData;
			this.targetData = targetData;
			this.corrSign = corrSign;
			this.mediators = mediators;
		}

		public String source;
		public String target;
		public SignedType edgeType;
		public RPPAData sourceData;
		public RPPAData targetData;
		public int corrSign;
		public String mediators;

		public String getEdgeData()
		{
			return source + "\t" + edgeType.getTag() + "\t" + target + "\t" + mediators;
		}
	}
}
