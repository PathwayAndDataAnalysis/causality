package org.cbio.causality.model;

import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.StudentsT;
import org.cbio.causality.util.Summary;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class RPPAData implements Cloneable
{
	public String id;
	public double[] vals0;
	public double[] vals1;
	public Set<String> genes;
	public Set<String> sites;
	public SiteEffect effect;
	private ChangeDetector chDet;

	public RPPAData(String id, double[] vals0, double[] vals1, Set<String> genes, Set<String> sites)
	{
		this(id, vals0, genes, sites);
		this.vals1 = vals1;
	}

	public RPPAData(String id, double[] vals, Set<String> genes, Set<String> sites)
	{
		this.id = id;
		this.vals0 = vals;
		this.genes = genes;
		this.sites = sites;

		if (sites == null) return;

		for (String site : sites)
		{
			for (String gene : genes)
			{
				Integer eff = PhosphoSitePlus.getEffect(gene, site);
				if (eff != null)
				{
					effect = SiteEffect.getValue(eff);
					break;
				}
			}
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
		return sites != null && !sites.isEmpty();
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
		if (vals1 == null) throw new UnsupportedOperationException();
		return Math.log(Summary.mean(vals1) / Summary.mean(vals0)) / Math.log(2);
	}

	public double getSignificanceBasedVal()
	{
		if (vals1 == null) throw new UnsupportedOperationException();
		double pval = getTTestPval();
		double ch = Summary.mean(vals1) - Summary.mean(vals0);

		double sig = -Math.log(pval) / Math.log(2);
		if (ch < 0) sig *= -1;

		return sig;
	}

	public double getMeanVal()
	{
		if (vals1 != null) throw new UnsupportedOperationException();
		return Summary.mean(vals0);
	}

	public double getLog2MeanVal()
	{
		if (vals1 != null) throw new UnsupportedOperationException();
		return Math.log(Summary.mean(vals0)) / Math.log(2);
	}

	public void separateData(boolean[][] subset)
	{
		if (vals1 != null) throw new UnsupportedOperationException();
		if (vals0.length != subset[0].length || vals0.length != subset[1].length)
			throw new IllegalArgumentException("Sizes don't match: vals0.length = " +
				vals0.length + ", subset[0].length = " + subset[0].length);

		double[] v0 = getValSubset(subset[0]);
		this.vals1 = getValSubset(subset[1]);
		this.vals0 = v0;
	}

	private double[] getValSubset(boolean[] sub)
	{
		double[] v = new double[ArrayUtil.countValue(sub, true)];
		int k = 0;
		for (int i = 0; i < sub.length; i++)
		{
			if (sub[i]) v[k++] = vals0[i];
		}
		return v;
	}

	public double getTTestPval()
	{
		if (vals1 == null) throw new UnsupportedOperationException();
		return StudentsT.getPValOfMeanDifference(vals0, vals1);
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

	@Override
	public Object clone()
	{
		try
		{
			return super.clone();
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
}
