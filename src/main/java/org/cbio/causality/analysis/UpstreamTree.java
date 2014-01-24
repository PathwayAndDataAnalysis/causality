package org.cbio.causality.analysis;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class UpstreamTree
{
	private Graph trav;
	private Graph lastStep;

	public UpstreamTree(Graph trav, Graph lastStep)
	{
		this.trav = trav;
		this.lastStep = lastStep;
	}

	public UpstreamTree(Graph trav)
	{
		this.trav = trav;
	}

	public GeneBranch getTree(String to, Set<String> from, int limit)
	{
		GeneBranch result = new GeneBranch(to);
		Map<String, GeneBranch> map = new HashMap<String, GeneBranch>();

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
			GeneBranch node = new GeneBranch(up);
			result.branches.add(node);
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
							map.get(gene).branches.add(map.get(up));
						}

						continue;
					}

					dist.put(up, i);
					visited.add(up);
					queue.add(up);

					GeneBranch node = new GeneBranch(up);

					if (!map.containsKey(gene))
					{
						System.out.println();
					}

					map.get(gene).branches.add(node);
					map.put(up, node);
				}
			}
		}

		result.removeNonTargetLeaf(from);
		return result;
	}

}
