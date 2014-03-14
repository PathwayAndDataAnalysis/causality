package org.cbio.causality.analysis;

import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import org.w3c.dom.Document;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author Emek Demir
 * @author Ozgun Babur
 */
public class RadialInfluenceTree
{

	private List<List<Node>> layers;

	private int rmult;

	private int smult;

	private int radiusOffset;

	private int stroke;

	private boolean flowToCenter;

	int xcnt = 700;

	int ycnt = 700;

	private Rectangle2D clip;

	public RadialInfluenceTree(List<List<Node>> layers, int rmult, int smult, int radiusOffset,
		int stroke, boolean flowToCenter)
	{
		this.layers = layers;
		this.rmult = rmult;
		this.smult = smult;
		this.radiusOffset = radiusOffset;
		this.stroke = stroke;
		this.flowToCenter = flowToCenter;
	}

	public RadialInfluenceTree(GeneBranch tree, boolean flowToCenter)
	{
		this(new ArrayList<List<Node>>(), 200, 20, 200, 20, flowToCenter);

		Map<String, Node> nodeMap = new HashMap<String, Node>();

		// create levels

		List<List<GeneBranch>> treeLevels = tree.getLevels();

		for (List<GeneBranch> treeLevel : treeLevels)
		{
			List<Node> layer = new ArrayList<Node>();

			for (GeneBranch branch : treeLevel)
			{
				Node node = createOrFind(branch, nodeMap);
				layer.add(node);
			}
			layers.add(layer);
		}

		// connect

		for (List<GeneBranch> treeLevel : treeLevels)
		{
			for (GeneBranch branch : treeLevel)
			{
				Node node = nodeMap.get(branch.gene);

				for (GeneBranch outer : branch.branches)
				{
					Node oNode = nodeMap.get(outer.gene);
					if(!oNode.out.contains(node))
					{
						oNode.out.add(node);
						node.in.add(oNode);
					}
				}
			}
		}

		merge();
	}

	public static void write(GeneBranch tree, boolean flowToCenter, String filename)
	{
		RadialInfluenceTree rit = new RadialInfluenceTree(tree, flowToCenter);
		rit.layout();
		try
		{
			String svgNS = "http://www.w3.org/2000/svg";

			Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().
				getDOMImplementation().createDocument(svgNS, "svg", null);

			SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(doc);
			ctx.setEmbeddedFontsOn(true);
			ctx.setPrecision(3);
			SVGGraphics2D svgGraphics2d = new SVGGraphics2D(ctx, true);
			rit.draw(svgGraphics2d);
			FileWriter writer = new FileWriter(filename);
			svgGraphics2d.stream(writer);
			writer.flush();
		}
		catch (ParserConfigurationException e)
		{
			e.printStackTrace();
		}
		catch (SVGGraphics2DIOException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public RadialInfluenceTree(String file) throws IOException
	{
		this(new ArrayList<List<Node>>(), 200, 20, 200, 20, true);

		InputStream is = RadialInfluenceTree.class.getClassLoader().getResourceAsStream(file);
		if (is == null) is = new FileInputStream(file);

		BufferedReader reader = new BufferedReader(
				new InputStreamReader(is));

		String line = reader.readLine();
		Map<String, Node> nodeMap = new HashMap<String, Node>();
		while (line != null)
		{
			StringTokenizer tk = new StringTokenizer(line);
			Color bgrnd = new Color((int) (128 * (Math.random() + 1)), ((int) (128 * (Math.random() + 1))),
			                        (int) (128 * (Math.random() + 1)), 230);
			Node source = createOrFind(tk.nextToken(), nodeMap, bgrnd);
			String type = tk.nextToken();
			bgrnd = new Color((int) (128 * (Math.random() + 1)), ((int) (128 * (Math.random() + 1))),
			                  (int) (128 * (Math.random() + 1)), 230);
			Node target = createOrFind(tk.nextToken(), nodeMap, bgrnd);
			if (!source.out.contains(target))
			{
				source.out.add(target);
				target.in.add(source);
			}
			line = reader.readLine();
		}
		Node root = nodeMap.get(file.substring(file.indexOf("/") + 1, file.indexOf(".")));
		List<Node> layer = new ArrayList<Node>();
		layer.add(root);
		this.layers.add(layer);
		findNextLayer(layer, new HashSet<Node>());
		merge();
	}

	private void merge()
	{
		List<Node> layer;
		for (int i = layers.size() - 1; i >= 0; i--)
		{
			layer = layers.get(i);
			HashMap<Node, Node> merge = new HashMap<Node, Node>();
			for (Node node : layer)
			{
				if (!merge.containsValue(node)) for (Node other : layer)
				{
					if (node != other && !node.in.isEmpty() && !node.out.isEmpty() && node.in.equals(other.in) &&
					    node.out.equals(other.out))
					{
						merge.put(node, other);
					}
				}
			}
			for (Node node : merge.keySet())
			{
				Node target = merge.get(node);
				target.name.putAll(node.name);
				layer.remove(node);
				for (Node node1 : node.in)
				{
					node1.out.remove(node);
				}
				for (Node node1 : node.out)
				{
					node1.in.remove(node);
				}
			}
			for (Node node : layer)
			{
				if (node.stroke == 0)
				{
					node.stroke = (int) (Math.random() * 3);
					for (Node in : node.in)
					{
						node.stroke += in.stroke;
					}
				}
			}
		}
	}

	private void findNextLayer(List<Node> layer, Set<Node> placed)
	{
		ArrayList<Node> next = new ArrayList<Node>();
		for (Node prev : layer)
		{
			for (Node node : prev.in)
			{
				if (!placed.contains(node)) next.add(node);
				placed.add(node);
			}
		}
		if (!next.isEmpty())
		{
			layers.add(next);
			findNextLayer(next, placed);
		}
	}


	private Node createOrFind(String s, Map<String, Node> nodeMap, Color color)
	{
		Node node = nodeMap.get(s);
		if (node == null)
		{
			node = new Node(s, (float) (1+3*Math.random()), new ArrayList<Node>(), new ArrayList<Node>(), color);
			nodeMap.put(s, node);
		}
		return node;
	}

	private Node createOrFind(GeneBranch branch, Map<String, Node> nodeMap)
	{
		Node node = nodeMap.get(branch.gene);
		if (node == null)
		{
			node = new Node(branch.gene, (float) (branch.getThickness()), new ArrayList<Node>(),
				new ArrayList<Node>(), branch.getColor());
			nodeMap.put(branch.gene, node);
		}
		return node;
	}

	public void draw(Graphics2D g2)
	{
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
//		g2.setClip(clip);
		for (int i = 1; i < layers.size(); i++)
		{
			List<Node> layer = layers.get(i);

			int tr = rmult * i;
			int trd = rmult * (i - 1);

			Color bgrnd = new Color(10, 10, 10, 10);
			g2.setColor(bgrnd);
			g2.fillOval(xcnt - tr, ycnt - tr, 2 * rmult * i, 2 * rmult * i);
			g2.setColor(Color.black);
			for (Node node : layer)
			{
				for (Node out : node.out)
				{
					Stroke old = g2.getStroke();
					g2.setStroke(new BasicStroke(node.stroke));
					if (i > 1)
					{
						QuadCurve2D curve = new QuadCurve2D.Double(
							node.getX(tr), node.getY(tr),
							out.getX(tr - rmult / 2), out.getY(tr - rmult / 2),
							out.getX(trd), out.getY(trd));

						g2.draw(curve);
						drawArrowToCurve(g2, curve);

					} else
					{
						double y1 = node.getY(tr);
						double x1 = node.getX(tr);
						g2.draw(new Line2D.Double(x1, y1, xcnt, ycnt));
						drawArrow(g2, (x1 + xcnt) / 2, (y1 + ycnt) / 2, y1 - ycnt, x1 - xcnt,
							Math.signum(ycnt - y1));
					}
					g2.setStroke(old);
				}
			}

		}

		drawLabels(g2);
	}

	private void drawArrowToCurve(Graphics2D g2, QuadCurve2D curve)
	{
		double cntrx = (curve.getX1() + curve.getX2() + 2 * curve.getCtrlX()) / 4;
		double cntry = (curve.getY1() + curve.getY2() + 2 * curve.getCtrlY()) / 4;
		double dy = (curve.getY1() + curve.getCtrlY()) / 2 - cntry;
		double dx = (curve.getX1() + curve.getCtrlX()) / 2 - cntrx;
		double dir = Math.signum(curve.getY2() - curve.getY1());
		drawArrow(g2, cntrx, cntry, dy, dx, dir);
	}

	private void drawArrow(Graphics2D g2, double cntrx, double cntry, double dy, double dx, double dir)
	{
		double xleft;
		double yleft;
		double xright;
		double yright;
		double arrowsize = 14;

		if (dy == 0)
		{
			yleft = cntry - 0.5 * arrowsize * Math.signum(dx);
			yright = cntry + 0.5 * arrowsize * Math.signum(dx);
			xleft = cntrx + arrowsize * Math.signum(dx);
			xright = cntrx + arrowsize * Math.signum(dx);

		}
		else if (dx == 0)
		{
			xleft = cntrx - 0.5 * arrowsize * Math.signum(dy);
			xright = cntrx + 0.5 * arrowsize * Math.signum(dy);
			yleft = cntry + arrowsize * Math.signum(dy);
			yright = cntry + arrowsize * Math.signum(dy);
		} else
		{
			double slopep = dy / dx;
			double slope = -dx / dy;

			double c = cntry - slope * cntrx;
			double signsiz = arrowsize * Math.signum(slopep) * dir;

			double arrx = 0.5 * signsiz / Math.sqrt(slope * slope + 1);
			double arrxp = signsiz / Math.sqrt(slopep * slopep + 1);

			xleft = (cntrx - arrx - arrxp);
			yleft = ((cntrx - arrx) * slope) + c - (arrxp * slopep);
			xright = (cntrx + arrx - arrxp);
			yright = ((cntrx + arrx) * slope) + c - (arrxp * slopep);
		}
		g2.draw(new Line2D.Double(xleft, yleft, cntrx, cntry));

		g2.draw(new Line2D.Double(xright, yright, cntrx, cntry));
	}

	private void drawLabels(Graphics2D g2)
	{
		Font mid = new Font(g2.getFont().getFamily(), 0, g2.getFont().getSize() * 2);
		Font large = new Font(g2.getFont().getFamily(), 0, g2.getFont().getSize() * 3);

		for (int i = 0; i < layers.size(); i++)
		{
			g2.setFont(i == 0 ? large : mid);
			FontMetrics metrics = g2.getFontMetrics();
			int exampleNameSize = metrics.stringWidth("ABCB1");
			List<Node> layer = layers.get(i);

			int tr = rmult * i;

			for (Node node : layer)
			{

				int p = metrics.getHeight();
				double min = Math.sqrt(exampleNameSize * node.name.size() * p);
				double current = 0;
				int line = 0;
				int size = node.name.keySet().size();
				double[] x = new double[size];
				double[] y = new double[size];
				double[] w = new double[size];
				List<String> sorted = new ArrayList(node.name.keySet());
				Collections.sort(sorted);
				for (int j = 0; j < sorted.size(); j++)
				{
					String s = sorted.get(j);
					if (current+5 + w[j] > min)
					{
						line++;
						current = 0;
					}
					x[j] = current;
					y[j] = line * p;
					w[j] = metrics.stringWidth(s);
					current += 5 + w[j];

				}
				double xoff= node.getX(tr)-(min+20)/2;
				double yoff= node.getY(tr)-p*(line+1)/2;
				for (int j = 0; j < sorted.size(); j++)
				{
					String s = sorted.get(j);
					g2.setColor(node.name.get(s));
					g2.fill(new RoundRectangle2D.Double(xoff + x[j] - 3, (yoff + y[j]), w[j] + 6, p, 10, 10));
					g2.setColor(Color.black);
					g2.drawString(s, ((float) (xoff + x[j])), (float) (yoff+y[j]+0.8*p));

				}
			}
		}
	}


	public void layout()
	{
		uniform();
		for (int i = 0; i < 30; i++)
		{
			upanddown(1);
		}
		for (int i = 1; i < 30; i++)
		{
			upanddown(i);
		}
		List<Node> all = getAllNodes();
		setClip();

		//rotateToTop();
	}

	private void rotateToTop(List<Node> all)
	{

		double center = findCenter(all);
		double shift = 1.5 * Math.PI - center % Math.PI;
		for (List<Node> layer : layers)
		{
			for (Node node : layer)
			{
				node.setAng(node.ang + shift);
			}
		}

	}

	private void setClip()
	{
		double xmax=Double.MIN_VALUE,ymax=Double.MIN_VALUE;
		double xmin=Double.MAX_VALUE,ymin=Double.MAX_VALUE;
		for (int i = 1; i < layers.size(); i++)
		{

			List<Node> layer = layers.get(i);
			int tr = rmult * i;
			for (Node node : layer)
			{
				double x = node.getX(tr);
				double y = node.getY(tr);
				if(x >xmax) xmax= x;
				if(x <xmin) xmin= x;
				if(y >ymax) ymax= y;
				if(y <ymin) ymin= y;
			}
		}
		this.xcnt= (int) (xcnt-xmin);
		this.ycnt= (int) (ycnt-ymin);
		this.clip = new Rectangle2D.Double(xmin-200,ymin-200,400+xmax-xmin,400+ymax-ymin);
	}

	private List<Node> getAllNodes()
	{
		List<Node> all = new ArrayList<Node>();
		for (int i = 1; i < layers.size(); i++)
		{
			List<Node> layer = layers.get(i);
			for (Node node : layer)
			{
				all.add(node);
			}
		}
		return all;
	}

	private void upanddown(double force)
	{
		for (int i = 1; i < layers.size(); i++)
		{
			List<Node> current = layers.get(i);
			Collections.sort(current);
			if (i == 1)
			{
				space(current, 1);
			} else
			{
				for (Node node : current)
				{

					if (node.out.size() > 0)
					{
						node.moveTowards(findCenter(node.out), 1.0 / force);
					}
				}
				space(current, layers.indexOf(current));
			}
		}
		for (int i = layers.size() - 1; i > 0; i--)
		{
			List<Node> current = layers.get(i);
			Collections.sort(current);
			for (Node node : current)
			{
				if (node.in.size() > 0)
				{
					node.moveTowards(findCenter(node.in), 1.0 / force);
				}

			}
			space(current, layers.indexOf(current));
		}

	}

	private double findCenter(List<Node> out)
	{
		if (out.size() == 1)
		{
			return out.get(0).ang;
		}
		Collections.sort(out);
		double cntr =((2 * Math.PI+ out.get(0).ang + out.get(out.size() - 1).ang)/ 2)%(2*Math.PI);
		double max = 2 * Math.PI + out.get(0).ang - out.get(out.size() - 1).ang ;
		for (int i = 0; i < out.size() - 1; i++)
		{
			double arc = out.get(i + 1).ang - out.get(i).ang;
			if (arc > max)
			{
				max = arc;
				cntr = ((((out.get(i + 1).ang + out.get(i).ang)) / 2));
			}
		}

		return cntr>Math.PI?cntr-Math.PI:cntr+Math.PI;

	}

	private void space(List<Node> current, int i)
	{
		Collections.sort(current);
		double mingap = Math.PI / (Math.sqrt(current.size() * 4) * i);   //TODO:determine mingap based on node sizes
		shift(current, mingap, 0);
	}

	private void shift(List<Node> current, double mingap, int i)
	{
		if (current.size() == 1) return;
		if (current.size() == i + 1)
		{
			double gap = (2 * Math.PI + current.get(0).ang - current.get(i).ang) % (2 * Math.PI);
//			if (gap < 0)
//			{
//
//				current.get(i).setAng(current.get(0).ang);
//				current.get(0).setAng(current.get(0).ang - gap);
//				Collections.sort(current);
//			}
			if (gap < mingap)
			{
				current.get(0).setAng(current.get(0).ang + (mingap - gap) / 2); //not perfect but will do ..
				current.get(i).setAng(current.get(i).ang - (mingap - gap) / 2); //not perfect but will do ..

			}
		} else
		{
			double gap = current.get(i + 1).ang - current.get(i).ang;
			if (gap < 0)
			{
				current.get(i).setAng(current.get(i + 1).ang);
				current.get(i + 1).setAng(current.get(i + 1).ang - gap / 2);
				Collections.sort(current);
			} else if (gap < mingap)
			{
				current.get(i + 1).setAng(current.get(i + 1).ang + (mingap - gap) / 2);
				current.get(i).setAng(current.get(i).ang - (mingap - gap) / 2);
			}
			shift(current, mingap, i + 1);
		}
	}


	private void uniform()
	{
		for (List<Node> layer : layers)
		{
			uniform(layer);
		}
	}

	private void uniform(List<Node> layer)
	{
		double gap = 2 * Math.PI / layer.size();
		for (int i = 0; i < layer.size(); i++)
		{
			layer.get(i).setAng(i * gap);
		}
	}


	public class Node implements Comparable<Node>
	{
		private HashMap<String, Color> name;

		private double ang = 0;

		private float stroke;

		private List<Node> out;

		private List<Node> in;

		public Node(String name,  float stroke, List<Node> out, List<Node> in, Color color)
		{
			this.name = new HashMap<String, Color>();
			this.name.put(name, color);

			this.stroke = stroke;
			this.out = out;
			this.in = in;
		}

		public double getX(double r)
		{

			return xcnt + (r * Math.cos(ang));
		}

		public double getY(double r)
		{
			return ycnt +  (r * Math.sin(ang));
		}

		public int compareTo(Node o)
		{
			return new Double(this.ang).compareTo(o.ang);
		}

		private void moveTowards(double ang, double force)
		{
			double delta = (ang - this.ang);
			if (delta > 0 && delta < Math.PI)
			{
				this.setAng(this.ang + delta * force);
			} else if (delta > 0 && delta > Math.PI)
			{
				this.setAng(this.ang - (delta - Math.PI) * force);
			} else if (delta < 0 && delta > -Math.PI)
			{
				this.setAng(this.ang + delta * force);
			} else if (delta < 0 && delta < -Math.PI)
			{
				this.setAng(this.ang - (delta + Math.PI) * force);
			}
		}

		private void setAng(double ang)
		{

			this.ang = ((2 * Math.PI) + ang) % (2 * Math.PI);

		}

		@Override
		public String toString()
		{
			return name.toString();
		}
	}

	public static void main(String[] args) throws IOException
	{
		BranchDataProvider dp = new BranchDataProvider()
		{
			@Override
			public Color getColor(String gene, String root)
			{
				return Color.WHITE;
			}

			@Override
			public double getThickness(GeneBranch branch, String root)
			{
				return 2;
			}
		};
		String r = "root";
		GeneBranch tree = new GeneBranch(r, dp, r);
		int i = 1;
		int n = 3;
		String prefix = "";
		for (int j = 0; j < n; j++) tree.branches.add(new GeneBranch(prefix + (i++) , dp, r));
		for (GeneBranch br : tree.branches)
		{
			for (int j = 0; j < n; j++) br.branches.add(new GeneBranch(prefix + (i++) , dp, r));

			for (GeneBranch brr : br.branches)
			{
				for (int j = 0; j < n; j++) brr.branches.add(new GeneBranch(prefix + (i++) , dp, r));
			}
		}

		write(tree, true, "temp.svg");
	}

//	public static void main(String[] args) throws IOException
//	{
//
//		final RadialInfluenceTree rit = new RadialInfluenceTree("sifs/TMEM164.sif");
////		final RadialInfluenceTree rit = new RadialInfluenceTree("temp/SH3BP4.sif");
//		rit.layout();
//		SwingUtilities.invokeLater(new Runnable()
//		{
//			public void run()
//			{
//				createAndShowGUI(rit);
//			}
//		});
//	}

	private static void createAndShowGUI(final RadialInfluenceTree rit)
	{
		JFrame f = new JFrame("Radial Influence Tree demo");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.add(new JPanel()
		{

			public Dimension getPreferredSize()
			{
				return new Dimension(1000, 1000);
			}

			public void paintComponent(Graphics g)
			{

				super.paintComponent(g);

				rit.draw((Graphics2D) g);

			}

		});
		f.pack();
		f.setSize(2500, 2500);
		f.setVisible(true);
	}
}