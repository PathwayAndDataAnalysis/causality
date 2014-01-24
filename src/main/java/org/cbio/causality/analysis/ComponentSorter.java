package org.cbio.causality.analysis;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ComponentSorter
{
	private Graph graph;
	private List<Set<String>> components;

	public ComponentSorter(Set<String> genes, Graph graph)
	{
		this.graph = graph;
		this.components = new ArrayList<Set<String>>();
		for (String gene : genes)
		{
			add(gene);
		}
		sort();
	}

	public List<Set<String>> getComponents(int sizeThr)
	{
		List<Set<String>> result = new ArrayList<Set<String>>();
		for (Set<String> comp : components)
		{
			if (comp.size() >= sizeThr) result.add(comp);
		}
		return result;
	}

	private void add(String gene)
	{
		Set<String> neighbors = graph.getNeighbors(gene);

		List<Set<String>> matches = new ArrayList<Set<String>>();

		for (Set<String> comp : components)
		{
			if (intersects(comp, neighbors)) matches.add(comp);
		}

		if (matches.isEmpty())
		{
			Set<String> set = new HashSet<String>();
			set.add(gene);
			components.add(set);
		}
		else if (matches.size() == 1)
		{
			matches.get(0).add(gene);
		}
		else
		{
			Set<String> set = matches.get(0);
			set.add(gene);
			for (int i = 1; i < matches.size(); i++)
			{
				set.addAll(matches.get(i));
			}
			for (int i = matches.size() - 1; i > 0; i--)
			{
				components.remove(matches.get(i));
			}
		}
	}

	private void sort()
	{
		Collections.sort(components, new Comparator<Set<String>>()
		{
			@Override
			public int compare(Set<String> o1, Set<String> o2)
			{
				return new Integer(o2.size()).compareTo(o1.size());
			}
		});
	}

	private boolean intersects(Set<String> set1, Set<String> set2)
	{
		Set<String> temp = new HashSet<String>(set1);
		temp.retainAll(set2);
		return !temp.isEmpty();
	}
}
