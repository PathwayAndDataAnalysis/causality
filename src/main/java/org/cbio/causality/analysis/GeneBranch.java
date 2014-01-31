package org.cbio.causality.analysis;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
* @author Ozgun Babur
*/
public class GeneBranch
{
	public String gene;
	public String root;
	public Set<GeneBranch> branches;
	boolean selected = false;
	public Set<GeneBranch> equivalents;
	private BranchDataProvider data;

	GeneBranch(String gene, BranchDataProvider data, String root)
	{
		this.gene = gene;
		this.data = data;
		this.root = root;
		this.branches = new HashSet<GeneBranch>();
	}

	GeneBranch(String gene, BranchDataProvider data)
	{
		this(gene, data, gene);
	}

	public boolean isLeaf()
	{
		return this.branches.isEmpty();
	}

	public boolean isSelected()
	{
		return selected;
	}

	public boolean selectLeaves(Collection<String> leaves)
	{
		if (isLeaf())
		{
			return selected = leaves.contains(gene);
		}

		boolean modified = false;

		for (GeneBranch branch : branches)
		{
			modified = branch.selectLeaves(leaves) || modified;
		}

		return this.selected = modified;
	}

	public void deselect()
	{
		this.selected = false;
		for (GeneBranch branch : branches)
		{
			branch.deselect();
		}
	}

	public boolean removeNonTargetLeaf(Collection<String> targets)
	{
		boolean modified = false;

		for (GeneBranch gwu : new HashSet<GeneBranch>(branches))
		{
			if (!gwu.branches.isEmpty())
			{
				boolean mod;
				do
				{
					mod = gwu.removeNonTargetLeaf(targets);
					modified = mod || modified;
				} while(mod);
			}
			if (gwu.branches.isEmpty() && !targets.contains(gwu.gene))
			{
				branches.remove(gwu);
				modified = true;
			}
		}

		return modified;
	}

	public void trimToMajorPaths(Collection<String> targets)
	{
		Map<GeneBranch, Set<String>> contentMap = new HashMap<GeneBranch, Set<String>>();

		for (GeneBranch u : branches)
		{
			Set<String> genes = u.getAllGenes();
			genes.retainAll(targets);
			contentMap.put(u, genes);
		}

		Set<GeneBranch> remove = new HashSet<GeneBranch>();

		for (GeneBranch g1 : contentMap.keySet())
		{
			for (GeneBranch g2 : contentMap.keySet())
			{
				if (g1 == g2) continue;

				if (contentMap.get(g1).containsAll(contentMap.get(g2)) &&
					contentMap.get(g1).size() > contentMap.get(g2).size())
				{
					remove.add(g2);
				}
			}
		}

		branches.removeAll(remove);

		for (GeneBranch u : branches)
		{
			u.trimToMajorPaths(targets);
		}
	}

	public boolean removeAbsentLeafConnections(Graph graph, boolean upstream)
	{
		boolean modified = false;

		for (GeneBranch branch : new HashSet<GeneBranch>(branches))
		{
			if (!branch.branches.isEmpty())
			{
				boolean mod;
				do
				{
					mod = branch.removeAbsentLeafConnections(graph, upstream);
					modified = mod || modified;
				} while(mod);
			}
			if (branch.branches.isEmpty() &&
				!(upstream ? graph.getUpstream(branch.gene) :
					         graph.getDownstream(branch.gene)).contains(gene))
			{
				branches.remove(branch);
				modified = true;
			}
		}

		return modified;
	}

	public Set<String> getAllGenes()
	{
		Set<String> set = new HashSet<String>();
		set.add(gene);
		for (GeneBranch u : branches)
		{
			set.addAll(u.getAllGenes());
		}
		return set;
	}

	public Set<String> getLeafGenes()
	{
		Set<String> set = new HashSet<String>();
		if (this.isLeaf()) set.add(gene);
		else
		{
			for (GeneBranch u : branches)
			{
				set.addAll(u.getLeafGenes());
			}
		}
		return set;
	}

	public Set<String> getAllGenesWithEq()
	{
		Set<String> set = getAllGenes();
		if (equivalents != null)
		{
			for (GeneBranch eq : equivalents)
			{
				set.add(eq.gene);
			}
		}
		return set;
	}

	private Map<GeneBranch, Set<GeneBranch>> getParentsMapping(Map<GeneBranch, Set<GeneBranch>> map)
	{
		if (!isLeaf())
		{
			if (map == null) map = new HashMap<GeneBranch, Set<GeneBranch>>();
			for (GeneBranch branch : branches)
			{
				if (!map.containsKey(branch)) map.put(branch, new HashSet<GeneBranch>());
				map.get(branch).add(this);
				branch.getParentsMapping(map);
			}
		}
		return map;
	}

	public void identifyEquivalents(boolean selectedOnly)
	{
		identifyEquivalents(null, selectedOnly);
	}

	private void identifyEquivalents(Map<GeneBranch, Set<GeneBranch>> parents, boolean selectedOnly)
	{
		if (isLeaf()) return;
		if (parents == null) parents = getParentsMapping(null);

		for (GeneBranch b1 : branches)
		{
			if (!b1.selected) continue;

			for (GeneBranch b2 : branches)
			{
				if (!b2.selected) continue;

				if (b2 != b1)
				{
					if (b1.branches.size() == b2.branches.size() &&
						b1.branches.containsAll(b2.branches))
					{
						Set<GeneBranch> p1 = parents.get(b1);
						Set<GeneBranch> p2 = parents.get(b2);

						if (p1.size() == p2.size() && p1.containsAll(p2))
						{
							if (b1.equivalents == null) b1.equivalents = new HashSet<GeneBranch>();
							if (b2.equivalents == null) b2.equivalents = new HashSet<GeneBranch>();

							b1.equivalents.add(b2);
							b2.equivalents.add(b1);
						}
					}
				}
			}
			b1.identifyEquivalents(parents, selectedOnly);
		}
	}

	public List<List<GeneBranch>> getLevels()
	{
		List<List<GeneBranch>> levels = new ArrayList<List<GeneBranch>>();

		List<GeneBranch> level = new ArrayList<GeneBranch>(1);
		level.add(this);

		Set<GeneBranch> visited = new HashSet<GeneBranch>();

		do
		{
			levels.add(level);

			// check edge sanity
			for (GeneBranch b : level) assert !visited.contains(b);
			visited.addAll(level);

			List<GeneBranch> next = new ArrayList<GeneBranch>();

			for (GeneBranch br : level)
			{
				next.addAll(br.branches);
			}

			level = next;
		}
		while (!level.isEmpty());

		return levels;
	}

	public Color getColor()
	{
		return data.getColor(gene, root);
	}

	public double getThickness()
	{
		return data.getThickness(this, root);
	}

	public GeneBranch copy(boolean selectedOnly)
	{
		if (selectedOnly && !this.isSelected()) return null;

		GeneBranch copy = new GeneBranch(gene, data, root);

		for (GeneBranch branch : branches)
		{
			GeneBranch bc = branch.copy(selectedOnly);
			if (bc != null) copy.branches.add(bc);
		}
		return copy;
	}

	@Override
	public String toString()
	{
		return gene + "\t" + branches.size();
	}
}
