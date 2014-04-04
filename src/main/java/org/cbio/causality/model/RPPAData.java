package org.cbio.causality.model;

import org.cbio.causality.network.PhosphoSitePlus;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class RPPAData
{
	public String id;
	public double[] vals;
	public Set<String> genes;
	public Set<String> sites;
	public SiteEffect effect;

	public RPPAData(String id, double[] vals, Set<String> genes, Set<String> sites)
	{
		this.id = id;
		this.vals = vals;
		this.genes = genes;
		this.sites = sites;

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

		if (effect == null && !sites.isEmpty())
		{
			for (String site : sites)
			{
				for (String gene : genes)
				{
					Integer eff = PhosphoSitePlus.getClosestEffect(gene, site);
					if (eff != null)
					{
						effect = SiteEffect.getValue(eff);
						break;
					}
				}
			}
		}
	}

	public boolean isPhospho()
	{
		return !sites.isEmpty();
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
		return obj instanceof RPPAData && this.genes.equals(((RPPAData) obj).genes);
	}

	@Override
	public int hashCode()
	{
		int code = 0;
		for (String gene : genes)
		{
			code += gene.hashCode();
		}
		return code;
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
}
