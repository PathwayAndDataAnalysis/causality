package org.cbio.causality.analysis;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class UpstreamTree
{
	private Traverse trav;
	private Traverse lastStep;

	public UpstreamTree(Traverse trav, Traverse lastStep)
	{
		this.trav = trav;
		this.lastStep = lastStep;
	}

	public UpstreamTree(Traverse trav)
	{
		this.trav = trav;
	}

	public GeneWithUp getTree(String to, Set<String> from, int limit)
	{
		GeneWithUp result = new GeneWithUp(to);
		Map<String, GeneWithUp> map = new HashMap<String, GeneWithUp>();

		// First step

		Map<String, Integer> dist = new HashMap<String, Integer>();
		Set<String> visited = new HashSet<String>();
		visited.add(to);
		dist.put(to, 0);
		LinkedList<String> queue = new LinkedList<String>();

		for (String up : lastStep != null ? lastStep.getUpstream(to) : trav.getUpstream(to))
		{
			dist.put(up, 1);
			visited.add(up);
			queue.add(up);
			GeneWithUp node = new GeneWithUp(up);
			result.up.add(node);
			map.put(up, node);
		}

		// Other steps

		for (int i = 2; i <= limit; i++)
		{
			if (queue.isEmpty()) break;

			List<String> temp = new ArrayList<String>(queue);
			queue.clear();

			for (String gene : temp)
			{
				assert visited.contains(gene);
				assert map.containsKey(gene);

				for (String up : trav.getUpstream(gene))
				{
					if (visited.contains(up))
					{
						if (dist.get(up) == dist.get(gene) + 1)
						{
							map.get(gene).up.add(map.get(up));
						}

						continue;
					}

					dist.put(up, i);
					visited.add(up);
					queue.add(up);

					GeneWithUp node = new GeneWithUp(up);

					if (!map.containsKey(gene))
					{
						System.out.println();
					}

					map.get(gene).up.add(node);
					map.put(up, node);
				}
			}
		}

		result.removeNonTargetLeaf(from);
		return result;
	}

	public class GeneWithUp
	{
		public String gene;
		public Set<GeneWithUp> up;

		GeneWithUp(String gene)
		{
			this.gene = gene;
			this.up = new HashSet<GeneWithUp>();
		}

		public boolean removeNonTargetLeaf(Set<String> targets)
		{
			boolean modified = false;

			for (GeneWithUp gwu : new HashSet<GeneWithUp>(up))
			{
				if (!gwu.up.isEmpty())
				{
					boolean mod;
					do
					{
						mod = gwu.removeNonTargetLeaf(targets);
						modified = mod || modified;
					} while(mod);
				}
				if (gwu.up.isEmpty() && !targets.contains(gwu.gene))
				{
					up.remove(gwu);
					modified = true;
				}
			}

			return modified;
		}

		public void trimToMajorPaths(Set<String> targets)
		{
			Map<GeneWithUp, Set<String>> contentMap = new HashMap<GeneWithUp, Set<String>>();

			for (GeneWithUp u : up)
			{
				Set<String> genes = u.getAllGenes();
				genes.retainAll(targets);
				contentMap.put(u, genes);
			}

			Set<GeneWithUp> remove = new HashSet<GeneWithUp>();

			for (GeneWithUp g1 : contentMap.keySet())
			{
				for (GeneWithUp g2 : contentMap.keySet())
				{
					if (g1 == g2) continue;

					if (contentMap.get(g1).containsAll(contentMap.get(g2)) &&
						contentMap.get(g1).size() > contentMap.get(g2).size())
					{
						remove.add(g2);
					}
				}
			}

			up.removeAll(remove);

			for (GeneWithUp u : up)
			{
				u.trimToMajorPaths(targets);
			}
		}

		public Set<String> getAllGenes()
		{
			Set<String> set = new HashSet<String>();
			set.add(gene);
			for (GeneWithUp u : up)
			{
				set.addAll(u.getAllGenes());
			}
			return set;
		}

		@Override
		public String toString()
		{
			return gene + "\t" + up.size();
		}
	}
}
