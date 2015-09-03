package org.cbio.causality.rppa;

import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.StudentsT;
import org.cbio.causality.util.Summary;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RPPAData implements Cloneable
{
	public String id;
	public double[][] vals;
	public List<String> genes;
	public Map<String, List<String>> sites;
	public SiteEffect effect;
	private ChangeDetector chDet;
	public String[][] header;
	Type type;

	public RPPAData(String id, double[][] vals, List<String> genes, Map<String, List<String>> sites)
	{
		this.id = id;
		this.vals = vals;
		this.genes = genes;
		this.sites = sites;

		if (sites != null && !sites.isEmpty()) type = Type.SITE_SPECIFIC;
		else type = Type.TOTAL_PROTEIN;

		if (sites == null) return;

		for (String gene : sites.keySet())
		{
			for (String site : sites.get(gene))
			{
				Integer eff = PhosphoSitePlus.getEffect(gene, site);
				if (eff != null)
				{
					effect = SiteEffect.getValue(eff);
					break;
				}
			}
			if (effect != null) break;
		}

//		if (effect == null && !sites.isEmpty())
//		{
//			for (String site : sites)
//			{
//				for (String gene : genes)
//				{
//					Integer eff = PhosphoSitePlus.getClosestEffect(gene, site);
//					if (eff != null)
//					{
//						effect = SiteEffect.getValue(eff);
//						break;
//					}
//				}
//			}
//		}
	}

	public void setChDet(ChangeDetector chDet)
	{
		this.chDet = chDet;
	}

	public boolean isPhospho()
	{
		return type == Type.SITE_SPECIFIC;
	}

	public boolean isTotalProt()
	{
		return type == Type.TOTAL_PROTEIN;
	}

	public boolean isActivity()
	{
		return type == Type.ACTIVITY;
	}

	public int getSelfEffect()
	{
		if (!isPhospho()) return 1;
		else if (effect == null) return 0;
		else return effect.getVal();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof RPPAData && id.equals(((RPPAData) obj).id);
	}

	@Override
	public int hashCode()
	{
		return id.hashCode();
	}

	public enum SiteEffect
	{
		ACTIVATING(1),
		INHIBITING(-1),
		COMPLEX(0);

		int val;

		private SiteEffect(int val)
		{
			this.val = val;
		}

		public int getVal()
		{
			return val;
		}

		public static SiteEffect getValue(int x)
		{
			for (SiteEffect effect : values())
			{
				if (effect.val == x) return effect;
			}
			return null;
		}
	}

	public double getLog2Ratio()
	{
		if (vals[1] == null) throw new UnsupportedOperationException();
		return Math.log(Summary.mean(vals[1]) / Summary.mean(vals[0])) / Math.log(2);
	}

	public double getSignificanceBasedVal()
	{
		if (vals[1] == null) throw new UnsupportedOperationException();
		double pval = getTTestPval();
		double ch = Summary.mean(vals[1]) - Summary.mean(vals[0]);

		double sig = -Math.log(pval) / Math.log(2);
		if (ch < 0) sig *= -1;

		return sig;
	}

	public double getMeanVal()
	{
		if (vals.length > 1 && vals[1] != null) throw new UnsupportedOperationException();
		return Summary.mean(vals[0]);
	}

	public double getLog2MeanVal()
	{
		if (vals[1] != null) throw new UnsupportedOperationException();
		return Math.log(Summary.mean(vals[0])) / Math.log(2);
	}

	public double getDifOfMeans()
	{
		if (vals[1] == null) throw new UnsupportedOperationException();
		return Summary.mean(vals[1]) - Summary.mean(vals[0]);
	}

	public void separateData(boolean[][] subset)
	{
		if (vals.length > 1) throw new UnsupportedOperationException();

		if (vals[0].length != subset[0].length || vals[0].length != subset[1].length)
			throw new IllegalArgumentException("Sizes don't match: vals[0].length = " +
				vals[0].length + ", subset[0].length = " + subset[0].length);

		for (int i = 1; i < subset.length; i++)
		{
			if (subset[0].length != subset[i].length)
				throw new IllegalArgumentException("Subset array sizes don't match: subset[0].length = " +
					subset[0].length + ", subset[" + i + "].length = " + subset[i].length);
		}


		double[] all = vals[0];
		vals = new double[subset.length][];
		for (int i = 0; i < vals.length; i++)
		{
			vals[i] = getValSubset(all, subset[i]);
		}
	}

	private double[] getValSubset(double[] all, boolean[] sub)
	{
		double[] v = new double[ArrayUtil.countValue(sub, true)];
		int k = 0;
		for (int i = 0; i < sub.length; i++)
		{
			if (sub[i]) v[k++] = all[i];
		}
		return v;
	}

	public double getTTestPval()
	{
		if (vals[1] == null) throw new UnsupportedOperationException();
		return StudentsT.getPValOfMeanDifference(vals[0], vals[1]);
	}

	public int getChangeSign()
	{
		if (chDet == null) throw new UnsupportedOperationException("Please set the change " +
			"detector (chDet) before calling this method.");

		return chDet.getChangeSign(this);
	}

	public double getChangeValue()
	{
		if (chDet == null) throw new UnsupportedOperationException("Please set the change " +
			"detector (chDet) before calling this method.");

		return chDet.getChangeValue(this);
	}

	public int getActvityChangeSign()
	{
		int eff = !isPhospho() ? 1 : effect == null ? 0 : effect.getVal();
		return getChangeSign() * eff;
	}

	public void makeActivityNode(boolean isActivated)
	{
		type = Type.ACTIVITY;
		vals = new double[][]{new double[1]};
		vals[0][0] = isActivated ? 1 : -1;
		setChDet(new DefaultActivityDet());
		effect = null;
		sites = null;
	}

	@Override
	public Object clone()
	{
		try
		{
			RPPAData clone = (RPPAData) super.clone();

//			clone.vals = vals.clone();

			return clone;
		}
		catch (CloneNotSupportedException e)
		{
			throw new RuntimeException(e);
		}
	}

	@Override
	public String toString()
	{
		return id;
	}

	public interface ChangeDetector
	{
		public int getChangeSign(RPPAData data);
		public double getChangeValue(RPPAData data);
	}

	public class DefaultActivityDet implements RPPAData.ChangeDetector
	{

		@Override
		public int getChangeSign(RPPAData data)
		{
			return (int) vals[0][0];
		}

		@Override
		public double getChangeValue(RPPAData data)
		{
			return vals[0][0];
		}
	}

	public static abstract class ChangeAdapter implements RPPAData.ChangeDetector
	{
		protected double threshold;

		@Override
		public int getChangeSign(RPPAData data)
		{
			double val = getChangeValue(data);
			if (val >= threshold) return 1;
			if (val <= -threshold) return -1;
			return 0;
		}

		@Override
		public double getChangeValue(RPPAData data)
		{
			return data.getMeanVal();
		}

		public void setThreshold(double threshold)
		{
			this.threshold = threshold;
		}
	}

	public static class TTestDetector extends ChangeAdapter
	{
		@Override
		public int getChangeSign(RPPAData data)
		{
			double pval = data.getTTestPval();
			if (pval > threshold) return 0;
			if (getChangeValue(data) > 0) return 1;
			return -1;
		}

		@Override
		public double getChangeValue(RPPAData data)
		{
			return data.getSignificanceBasedVal();
		}
	}

	public static void write(Collection<RPPAData> datas, String filename)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

			writer.write("ID\tSymbol\tSite\tEffect");

			Iterator<RPPAData> iter = datas.iterator();
			RPPAData sample = iter.next();
			while (sample.isActivity()) sample = iter.next();

			for (int j = 0; j < sample.vals.length; j++)
			{
				String[] header = sample.header == null ? null : sample.header[j];
				for (int i = 0; i < sample.vals[j].length; i++)
				{
					writer.write(header == null ? "\tv" + j + "-" + i : "\t" + header[i]);
				}
			}

			for (RPPAData data : datas)
			{
				if (data.isActivity()) continue;

				writer.write("\n" + data.id);

				String items = new ArrayList<String>(data.genes).toString().
					replace(",", "").replace("[", "").replace("]", "");
				writer.write("\t" + items);

				if (data.sites != null)
				{
					List<String> siteList = new ArrayList<String>();
					for (String gene : data.sites.keySet())
					{
						List<String> ss = data.sites.get(gene);
						String s = ss.toString().replace("[", "").replace("]", "").replace(", ", "|");
						siteList.add(s);
					}

					items = siteList.toString().replace(",", "").replace("[", "").replace("]", "");
					writer.write("\t" + items);

					writer.write("\t" + (data.effect == null ? "" : data.effect == SiteEffect.COMPLEX ?
						"c" : data.effect == SiteEffect.ACTIVATING ? "a" : "i"));
				}
				else writer.write("\t\t");

				for (double[] vals : data.vals)
				{
					for (double v : vals)
					{
						writer.write("\t" + v);
					}
				}
			}

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public enum Type
	{
		TOTAL_PROTEIN,
		SITE_SPECIFIC,
		ACTIVITY,
		EXPRESSION
	}

	public static void shuffleValues(Collection<RPPAData> datas)
	{
		List<RPPAData> list = new ArrayList<RPPAData>();

		for (RPPAData data : datas)
		{
			if (!data.isActivity()) list.add(data);
		}

		List<Integer> index = new ArrayList<Integer>();
		for (int i = 0; i < datas.size(); i++)
		{
			index.add(i);
		}
		Collections.shuffle(index);

		RPPAData d0 = list.get(index.get(0));
		double[][] v = d0.vals;

		for (int i = 1; i < index.size(); i++)
		{
			RPPAData dd = list.get(index.get(i));
			double[][] temp = dd.vals;
			dd.vals = v;
			v = temp;
		}

		d0.vals = v;
	}

	public static List<RPPAData> copy(List<RPPAData> orig)
	{
		List<RPPAData> list = new ArrayList<RPPAData>();
		for (RPPAData data : orig)
		{
			list.add((RPPAData) data.clone());
		}
		return list;
	}
}
