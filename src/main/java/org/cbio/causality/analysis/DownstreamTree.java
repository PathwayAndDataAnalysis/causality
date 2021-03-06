package org.cbio.causality.analysis;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class DownstreamTree
{
	private Graph trav;
	private Graph lastStep;
	private BranchDataProvider data;

	public DownstreamTree(Graph trav, Graph lastStep, BranchDataProvider data)
	{
		this.trav = trav;
		this.lastStep = lastStep;
		this.data = data;
	}

	public DownstreamTree(Graph trav, Graph lastStep)
	{
		this(trav, lastStep, null);
	}

	public DownstreamTree(Graph trav)
	{
		this(trav, null, null);
	}

	public GeneBranch getTree(String from, Collection<String> to, int limit)
	{
		GeneBranch result = new GeneBranch(from, data);
		Map<String, GeneBranch> map = new HashMap<String, GeneBranch>();
		map.put(from, result);

		Map<String, Integer> dist = new HashMap<String, Integer>();
		Set<String> visited = new HashSet<String>();
		Set<String> nonLeafVisited = new HashSet<String>();
		visited.add(from);
		nonLeafVisited.add(from);
		dist.put(from, 0);
		LinkedList<String> queue = new LinkedList<String>();
		queue.add(from);

		if (lastStep != null) expandLastStep(map, dist, visited, from, from);

		for (int i = 0; i < (lastStep == null ? limit : limit - 1); i++)
		{
			if (queue.isEmpty()) break;

			List<String> temp = new ArrayList<String>(queue);
			queue.clear();

			for (String gene : temp)
			{
				assert visited.contains(gene);
				assert nonLeafVisited.contains(gene);
				assert map.containsKey(gene);

				for (String down : trav.getDownstream(gene))
				{
					if (visited.contains(down))
					{
						if (dist.get(down) == dist.get(gene) + 1)
						{
							map.get(gene).branches.add(map.get(down));

							if (!nonLeafVisited.contains(down))
							{
								queue.add(down);
								nonLeafVisited.add(down);
							}
						}

						continue;
					}

					dist.put(down, dist.get(gene) + 1);
					queue.add(down);
					visited.add(down);
					nonLeafVisited.add(down);

					GeneBranch node = new GeneBranch(down, data, from);

					map.get(gene).branches.add(node);
					map.put(down, node);
				}
			}

			if (lastStep != null)
			{
				for (String gene : queue)
				{
					expandLastStep(map, dist, visited, gene, from);
				}
			}
		}

		result.removeNonTargetLeaf(to);
		if (lastStep != null) result.removeAbsentLeafConnections(lastStep, true);

		return result;
	}

	private void expandLastStep(Map<String, GeneBranch> map, Map<String, Integer> dist,
		Set<String> visited, String gene, String root)
	{
		for (String down : lastStep.getDownstream(gene))
		{
			if (visited.contains(down))
			{
				if (dist.get(down) == dist.get(gene) + 1)
				{
					GeneBranch node = map.get(down);
					map.get(gene).branches.add(node);
				}
			}
			else
			{
				dist.put(down, dist.get(gene) + 1);
				visited.add(down);
				GeneBranch node = new GeneBranch(down, data, root);
				map.get(gene).branches.add(node);
				map.put(down, node);
			}
		}
	}
}
