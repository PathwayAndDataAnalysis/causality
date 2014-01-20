package org.cbio.causality.analysis;

import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.io.*;
import java.util.*;

/**
 * This is a simple graph, built using external maps. All nodes are identified with a unique String.
 * Relations can be either directed or undirected, and can be mixed.
 *
 * @author Ozgun Babur
 */
public class Graph
{
	protected Map<String, Set<String>> dwMap;
	protected Map<String, Set<String>> upMap;
	protected Map<String, Set<String>> ppMap;

	private boolean allowSelfEdges = false;

	public Graph()
	{
		dwMap = new HashMap<String, Set<String>>();
		upMap = new HashMap<String, Set<String>>();
		ppMap = new HashMap<String, Set<String>>();
	}

	public boolean load(String filename, Set<String> ppiTypes, Set<String> signalTypes)
	{
		try
		{
			return load(new FileInputStream(filename), ppiTypes, signalTypes);
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			return false;
		}
	}

	public boolean load(InputStream is, Set<String> undirectedTypes, Set<String> directedTypes)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new InputStreamReader(is));

			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				String[] token = line.split("\t");

				if (token.length < 3) continue;

				if (undirectedTypes.contains(token[1]))
				{
					putRelation(token[0], token[2], false);
				}
				else if (directedTypes.contains(token[1]))
				{
					putRelation(token[0], token[2], true);
				}
			}

			reader.close();
		}
		catch (IOException e) { e.printStackTrace(); return false; } return true;
	}

	public void load(Collection<SIFInteraction> sifs, SIFType... typeArray)
	{
		Set<SIFType> types = typeArray.length == 0 ? null :
			new HashSet<SIFType>(Arrays.asList(typeArray));

		for (SIFInteraction sif : sifs)
		{
			if (types != null && !types.contains(sif.type)) continue;

			putRelation(sif.sourceID, sif.targetID, sif.type.isDirected());
		}
	}

	public void clear()
	{
		upMap.clear();
		dwMap.clear();
		ppMap.clear();
	}

	public void putRelation(String source, String target, boolean directed)
	{
		if (!allowSelfEdges && source.equals(target)) return;

		if (directed)
		{
			if (!upMap.containsKey(target)) upMap.put(target, new HashSet<String>());
			if (!dwMap.containsKey(source)) dwMap.put(source, new HashSet<String>());
			upMap.get(target).add(source);
			dwMap.get(source).add(target);
		}
		else
		{
			if (!ppMap.containsKey(source)) ppMap.put(source, new HashSet<String>());
			if (!ppMap.containsKey(target)) ppMap.put(target, new HashSet<String>());
			ppMap.get(source).add(target);
			ppMap.get(target).add(source);
		}
	}

	public Set<String> goBFS(String seed, boolean downstream)
	{
		return goBFS(seed, downstream ? dwMap : upMap);
	}

	public Set<String> goBFS(Set<String> seed, Set<String> visited, boolean downstream)
	{
		return goBFS(seed, visited, downstream ? dwMap : upMap);
	}

	public Set<String> goBFS(Set<String> seed, Set<String> visited)
	{
		return goBFS(seed, visited, ppMap);
	}

	private Set<String> goBFS(String seed, Map<String, Set<String>> map)
	{
		return goBFS(Collections.singleton(seed), null, map);
	}

	private Set<String> goBFS(Set<String> seed, Set<String> visited, Map<String, Set<String>> map)
	{
		Set<String> neigh = new HashSet<String>();
		for (String s : seed)
		{
			if (map.containsKey(s))
			{
				for (String n : map.get(s))
				{
					if (visited == null || !visited.contains(n))
					{
						neigh.add(n);
					}
				}
			}
		}
		return neigh;
	}

	public Set<String> getUpstream(Collection<String> genes)
	{
		Set<String> up = new HashSet<String>();
		for (String gene : genes)
		{
			up.addAll(getUpstream(gene));
		}
		return up;
	}

	public Set<String> getUpstream(String gene)
	{
		if (upMap.containsKey(gene)) return upMap.get(gene);
		else return Collections.emptySet();
	}

	public Set<String> getUpstream(String gene, int depth)
	{
		return getUpstream(Collections.singleton(gene), depth);
	}

	public Set<String> getUpstream(Set<String> genes, int depth)
	{
		if (depth < 1) throw new IllegalArgumentException(
			"Depth has to be positive and non-zero. Depth: " + depth);

		Set<String> newUp = new HashSet<String>(genes);
		Set<String> up = new HashSet<String>();

		for (int i = 0; i < depth; i++)
		{
			newUp = getUpstream(newUp);
			newUp.removeAll(up);
			up.addAll(newUp);
		}

		return up;
	}

	public Set<String> getDownstream(String gene)
	{
		if (dwMap.containsKey(gene)) return dwMap.get(gene);
		else return Collections.emptySet();
	}

	public Set<String> getNeighbors(String gene)
	{
		Set<String> n = new HashSet<String>();
		if (ppMap.get(gene) != null) n.addAll(ppMap.get(gene));
		if (upMap.get(gene) != null) n.addAll(upMap.get(gene));
		if (dwMap.get(gene) != null) n.addAll(dwMap.get(gene));
		return n;
	}
	
	public Set<String> getNeighbors(Set<String> genes)
	{
		Set<String> n = new HashSet<String>();
		for (String gene : genes)
		{
			n.addAll(getNeighbors(gene));
		}
		return n;
	}

	public int getDegree(String gene)
	{
		return getNeighbors(gene).size();
	}
	
	public Set<String> getPathElements(String from, Set<String> to, int limit)
	{
		Set<String> result = new HashSet<String>();
		getPathElements(from, to, limit, 0, result);
		return result;
	}

	private void getPathElements(String from, Set<String> to, int limit, int i, Set<String> result)
	{
		Set<String> set = Collections.singleton(from);
		Set<String> neigh = goBFS(set, set, true);
		for (String n : neigh)
		{
			if (to.contains(n)) result.add(n);
			else if (i < limit)
			{
				int prevSize = result.size();
				getPathElements(n, to, limit, i+1, result);
				if (result.size() > prevSize) result.add(n);
			}
		}
	}
	
	public List<CommPoint> getCommonDownstream(Set<String> seed, int limit)
	{
		Map<String, Set<String>> reachMap = new HashMap<String, Set<String>>();
		Map<String, Set<String>> breadthMap = new HashMap<String, Set<String>>();
		Map<String, Set<String>> visitedMap = new HashMap<String, Set<String>>();

		Set<CommPoint> points = new HashSet<CommPoint>();
		
		for (String s : seed)
		{
			reachMap.put(s, new HashSet<String>(Arrays.asList(s)));
			breadthMap.put(s, new HashSet<String>(Arrays.asList(s)));
			visitedMap.put(s, new HashSet<String>(Arrays.asList(s)));
		}

		for (int i = 1; i < limit; i++)
		{
			for (String s : seed)
			{
				Set<String> neigh = goBFS(breadthMap.get(s), visitedMap.get(s), true);
				for (String n : neigh)
				{
					if (!reachMap.containsKey(n))
						reachMap.put(n, new HashSet<String>(Arrays.asList(s)));
					else reachMap.get(n).add(s);
				}
				breadthMap.put(s, neigh);
				visitedMap.get(s).addAll(neigh);
			}

			for (String r : reachMap.keySet())
			{
				if (reachMap.get(r).size() > 1)
				{
					CommPoint p = new CommPoint(r, reachMap.get(r), i);
					if (!containsBetter(points, p)) points.add(p);
				}
			}
		}

		List<CommPoint> list = new ArrayList<CommPoint>(points);
		Collections.sort(list);
		return list;
	}
	
	private boolean containsBetter(Set<CommPoint> set, CommPoint p)
	{
		if (set.contains(p)) return true;
		for (CommPoint cp : set)
		{
			if (cp.dist < p.dist && cp.upstr.containsAll(p.upstr)) return true;
		}
		return false;
	}
	
	public Set<String> getSymbols()
	{
		Set<String> merge = new HashSet<String>(upMap.keySet());
		merge.addAll(dwMap.keySet());
		merge.addAll(ppMap.keySet());
		for (Set<String> vals : upMap.values()) merge.addAll(vals);
		for (Set<String> vals : dwMap.values()) merge.addAll(vals);
		for (Set<String> vals : ppMap.values()) merge.addAll(vals);
		return merge;
	}
	
	class CommPoint implements Comparable
	{
		String s;
		Set<String> upstr;
		int dist;

		CommPoint(String s, Set<String> upstr, int dist)
		{
			this.s = s;
			this.upstr = upstr;
			this.dist = dist;
		}

		@Override
		public boolean equals(Object o)
		{
			if (o instanceof CommPoint)
			{
				CommPoint p = (CommPoint) o;
				if (p.s.equals(s) && p.upstr.containsAll(upstr) && upstr.containsAll(p.upstr) && 
					p.dist == dist) return true;
			}
			return false;
		}

		@Override
		public int hashCode()
		{
			return s.hashCode();
		}

		@Override
		public int compareTo(Object o)
		{
			if (o instanceof CommPoint)
			{
				CommPoint p = (CommPoint) o;
				new Integer(p.upstr.size()).compareTo(upstr.size());
			}
			return 0;
		}
	}

	/**
	 * Gets the common downstream of length 1, but allows more length if path to downstream is also
	 * in the seed set.
	 */
	public Set<String> getLinkedCommonDownstream(Set<String> seed)
	{
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();

		for (String s : seed)
		{
			map.put(s, goBFS(Collections.singleton(s), null, true));
			map.get(s).add(s);
		}

		boolean loop = true;

		while(loop)
		{
			loop = false;

			for (String s1 : seed)
			{
				for (String s2 : seed)
				{
					if (s1.equals(s2)) continue;

					if (map.get(s2).contains(s1))
					{
						boolean changed = map.get(s2).addAll(map.get(s1));
						loop = changed || loop;
					}
				}
			}
		}

		Set<String> result = new HashSet<String>(map.values().iterator().next());

		for (Set<String> set : map.values())
		{
			result.retainAll(set);
		}

		return result;
	}

	public Map<Integer, Integer> getDegreeDistibution(boolean indegree)
	{
		Map<Integer, Integer> dist = new HashMap<Integer, Integer>();
		collectDegrees(dist, indegree ? upMap : dwMap);
		return dist;
	}

	public Map<Integer, Integer> getDegreeDistibution()
	{
		Map<Integer, Integer> dist = new HashMap<Integer, Integer>();
		collectDegrees(dist, upMap);
		collectDegrees(dist, dwMap);
		collectDegrees(dist, ppMap);
		return dist;
	}

	private void collectDegrees(Map<Integer, Integer> dist, Map<String, Set<String>> map)
	{
		for (Set<String> set : map.values())
		{
			int degree = set.size();

			if (!dist.containsKey(degree)) dist.put(degree, 1);
			else dist.put(degree, dist.get(degree) + 1);
		}
	}

	public static void main(String[] args)
	{
		Graph graph = new Graph();
		graph.load("/home/ozgun/Desktop/SIF.txt", new HashSet<String>(Arrays.asList("BINDS_TO")),
			new HashSet<String>(Arrays.asList("STATE_CHANGE", "TRANSCRIPTION", "DEGRADATION")));
		System.out.println("traverse.upMap.size() = " + graph.upMap.size());
		System.out.println("traverse.dwMap.size() = " + graph.dwMap.size());
		Set<String> merge = new HashSet<String>(graph.upMap.keySet());
		merge.addAll(graph.dwMap.keySet());
		System.out.println("merge.size() = " + merge.size());
		System.out.println("traverse.ppMap.size() = " + graph.ppMap.size());

		for (String n : graph.dwMap.get("EP300"))
		{
			System.out.println(n);
		}
	}
}
