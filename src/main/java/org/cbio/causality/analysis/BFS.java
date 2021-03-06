package org.cbio.causality.analysis;

import org.cbio.causality.model.Node;
import org.biopax.paxtools.query.algorithm.Direction;
import org.biopax.paxtools.query.model.Edge;
import org.biopax.paxtools.query.model.GraphObject;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

/**
 * Implements breadth-first search on the signed and directed graph.
 *
 * @author Ozgun Babur
 * @author Merve Cakir
 */
public class BFS
{
	/**
	 * Distance labels. Missing label interpreted as infinitive.
	 */
	protected Map<GraphObject, Integer> dist;

	/**
	 * Color labels. Missing color interpreted as white.
	 */
	protected Map<GraphObject, Integer> colors;

	/**
	 * BFS starts from source nodes. They get the label 0.
	 */
	protected Set<Node> sourceSet;

	/**
	 * BFS will not further traverse neighbors of any node in the stopSet.
	 */
	protected Set<Node> stopSet;

	/**
	 * Whether the direction is FORWARD, it is REVERSE otherwise.
	 */
	protected Direction direction;

	/**
	 * Stop distance.
	 */
	protected int limit;

	protected LinkedList<Node> queue;

	protected boolean startFromTranscription;
	
	public BFS(Set<Node> sourceSet, Set<Node> stopSet, Direction direction, int limit,
		boolean startFromTranscription)
	{
		if (direction == Direction.BOTHSTREAM)
			throw new IllegalArgumentException("BOTHSTREAM is not a valid parameter in BFS");

		if (startFromTranscription && direction == Direction.DOWNSTREAM)
			throw new IllegalArgumentException("Cannot start from transcription when going downstream");

		this.sourceSet = sourceSet;
		this.stopSet = stopSet;
		this.direction = direction;
		this.limit = limit;
		this.startFromTranscription = startFromTranscription;
	}

	/**
	 * Empty constructor for other possible uses.
	 */
	public BFS()
	{
	}

	public Map<GraphObject, Integer> run()
	{
		return runBFS();
	}

	public Map<GraphObject, Integer> runBFS()
	{
		initMaps();

		// Add all source nodes to the queue if traversal is needed

		if (limit > 0)
		{
			queue.addAll(sourceSet);
		}

		// Initialize dist and color of source set

		for (Node source : sourceSet)
		{
			setLabel(source, 0);
			setColor(source, GRAY);
			source.setPathSign(1);

			labelEquivRecursive(source, UPWARD, 0, true, false, 1);
			labelEquivRecursive(source, DOWNWARD, 0, true, false, 1);
		}

		// Process the queue

		while (!queue.isEmpty())
		{
			Node current = queue.remove(0);

			processNode(current);

			// Current node is processed
			setColor(current, BLACK);
		}

		return dist;
	}

	protected void initMaps()
	{
		// Initialize label, maps and queue

		dist = new HashMap<GraphObject, Integer>();
		colors = new HashMap<GraphObject, Integer>();
		queue = new LinkedList<Node>();
	}

	protected void processNode(Node current)
	{
		// Do not process the node if it is ubique

		if (current.isUbique())
		{
			setColor(current, BLACK);
			return;
		}

		// Process edges towards the direction

		for (Edge edge : direction == Direction.DOWNSTREAM ?
			current.getDownstream() : current.getUpstream())
		{
			assert edge != null;

			// Label the edge considering direction of traversal and type of current node

			if (direction == Direction.DOWNSTREAM || !current.isBreadthNode())
			{
				setLabel(edge, getLabel(current));
			}
			else
			{
				setLabel(edge, getLabel(current) + 1);
			}

			// Get the other end of the edge
			Node neigh = (Node) (direction == Direction.DOWNSTREAM ?
				edge.getTargetNode() : edge.getSourceNode());

			assert neigh != null;

			int pathSign = current.getPathSign() * edge.getSign();

			// Decide neighbor label according to the search direction and node type
			int dist = getLabel(edge);
			if (neigh.isBreadthNode() && direction == Direction.DOWNSTREAM) dist++;

			// Check if we need to stop traversing the neighbor, enqueue otherwise
			boolean further = (stopSet == null || !isEquivalentInTheSet(neigh, stopSet)) &&
				(!neigh.isBreadthNode() || dist < limit) && !neigh.isUbique() &&
				(!startFromTranscription || neigh.isBreadthNode() || dist > 1 || neigh.isTranscription()); // this line if for starting from transcriptions

			// Process the neighbor if not processed or not in queue

			if (getColor(neigh) == WHITE)
			{
				// Label the neighbor
				setLabel(neigh, dist);
				neigh.setPathSign(pathSign);

				if (further)
				{
					setColor(neigh, GRAY);

					// Enqueue the node according to its type

					if (neigh.isBreadthNode())
					{
						queue.addLast(neigh);
					}
					else
					{
						// Non-breadth nodes are added in front of the queue
						queue.addFirst(neigh);
					}
				}
				else
				{
					// If we do not want to traverse this neighbor, we paint it black
					setColor(neigh, BLACK);
				}

			}
			else if (getLabel(neigh) == dist)
			{
				if (neigh.getPathSign() != pathSign)
				{
					neigh.setPathSign(0);
				}
			}

			labelEquivRecursive(neigh, UPWARD, dist, further, !neigh.isBreadthNode(), pathSign);
			labelEquivRecursive(neigh, DOWNWARD, dist, further, !neigh.isBreadthNode(), pathSign);
		}
	}

	protected void labelEquivRecursive(Node node, boolean up, int dist,
		boolean enqueue, boolean head, int pathSign)
	{
		for (org.biopax.paxtools.query.model.Node equivalent : up ? 
			node.getUpperEquivalent() : node.getLowerEquivalent())
		{
			Node equiv = (Node) equivalent;
			
			if (getColor(equiv) == WHITE)
			{
				setLabel(equiv, dist);
				equiv.setPathSign(pathSign);
	
				if (enqueue)
				{
					setColor(equiv, GRAY);
	
					if (head) queue.addFirst(equiv);
					else queue.add(equiv);
				}
				else
				{
					setColor(equiv, BLACK);
				}
			}
			else if (getLabel(equiv) == dist)
			{
				if (equiv.getPathSign() != pathSign)
				{
					equiv.setPathSign(0);
				}				
			}

			labelEquivRecursive(equiv, up, dist, enqueue, head, pathSign);
		}
	}

	protected boolean isEquivalentInTheSet(Node node, Set<Node> set)
	{
		return set.contains(node) ||
			isEquivalentInTheSet(node, UPWARD, set) || isEquivalentInTheSet(node, DOWNWARD, set); 
	}

	protected boolean isEquivalentInTheSet(Node node, boolean direction, Set<Node> set)
	{
		for (org.biopax.paxtools.query.model.Node equivalent : direction == UPWARD ? 
			node.getUpperEquivalent() : node.getLowerEquivalent())
		{
			Node eq = (Node) equivalent;
			if (set.contains(eq)) return true;
			boolean isIn = isEquivalentInTheSet(eq, direction, set);
			if (isIn) return true;
		}
		return false;
	}

	protected int getColor(Node node)
	{
		if (!colors.containsKey(node))
		{
			// Absence of color is interpreted as white
			return WHITE;
		}
		else
		{
			return colors.get(node);
		}
	}

	protected void setColor(Node node, int color)
	{
		colors.put(node, color);
	}

	public int getLabel(GraphObject go)
	{
		if (!dist.containsKey(go))
		{
			// Absence of label is interpreted as infinite
			return Integer.MAX_VALUE-(limit*2);
		}
		else
		{
			return dist.get(go);
		}
	}

	protected void setLabel(GraphObject go, int label)
	{
		dist.put(go, label);
	}

	/**
	 * Towards parents.
	 */
	public static final boolean UPWARD = true;

	/**
	 * Towards children.
	 */
	public static final boolean DOWNWARD = false;

	/**
	 * Color white indicates the node is not processed.
	 */
	public static final int WHITE = 0;

	/**
	 * Color gray indicates that the node is in queue waiting to be processed.
	 */
	public static final int GRAY = 1;

	/**
	 * Color black indicates that the node was processed.
	 */
	public static final int BLACK = 2;
}
