package org.cbio.causality.alteredreaction;

/**
 * @author Ozgun Babur
 */
public class Interaction implements Comparable
{
	public String sym1;
	public String sym2;

	public Integer score;

	public Interaction(String sym1, String sym2)
	{
		if (sym1.compareTo(sym2) > 0)
		{
			this.sym2 = sym1;
			this.sym1 = sym2;
		}
		else
		{
			this.sym1 = sym1;
			this.sym2 = sym2;
		}

		assert sym1.compareTo(sym2) < 0;
		score = 0;
	}

	@Override
	public String toString()
	{
		return sym1 + "\t" + sym2;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof Interaction)
		{
			Interaction in = (Interaction) obj;
			return sym1.equals(in.sym1) && sym2.equals(in.sym2);
		}
		return false;
	}

	@Override
	public int hashCode()
	{
		return sym1.hashCode() + sym2.hashCode();
	}

	@Override
	public int compareTo(Object o)
	{
		if (o instanceof Interaction)
		{
			Interaction in = (Interaction) o;
			return in.score.compareTo(score);
		}
		return 0;
	}
}
