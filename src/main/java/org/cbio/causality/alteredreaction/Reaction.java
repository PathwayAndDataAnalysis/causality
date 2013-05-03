package org.cbio.causality.alteredreaction;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Reaction
{
	public String ID;
	public List<String> genes;
	Map<String, AlterationPack> alterations;
	Map<String, Integer> shares;

	public Change[] changes;
	Integer coverage;

	public Reaction()
	{
		genes = new ArrayList<String>();
		alterations = new HashMap<String, AlterationPack>();
	}

	public void fillShares()
	{
		shares = new HashMap<String, Integer>();
		for (String gene : genes)
		{
			AlterationPack pack = alterations.get(gene);
			int cnt = pack.getAlteredCount(Alteration.ANY);
			shares.put(gene, cnt);
		}
	}

	public void sortGenes()
	{
		if (shares == null) fillShares();

		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String gene1, String gene2)
			{
				int c = shares.get(gene2).compareTo(shares.get(gene1));
				if (c != 0) return c;
				return gene1.compareTo(gene2);
			}
		});
	}

	public void fillChanges()
	{
		changes = new Change[alterations.get(genes.get(0)).getSize()];
		int cnt = 0;

		for (int i = 0; i < changes.length; i++)
		{
			for (String gene : genes)
			{
				if (alterations.get(gene).getChange(Alteration.ANY, i).isAltered())
				{
					changes[i] = Change.ACTIVATING;
					cnt++;
					break;
				}
			}
			if (changes[i] == null) changes[i] = Change.NO_CHANGE;
		}

		coverage = cnt;
	}

	public boolean isSubsetOf(Reaction r)
	{
		return r.genes.containsAll(genes) && r.genes.size() > genes.size();
	}

	@Override
	public String toString()
	{
		String s = ID + "\t" + coverage;
		for (String gene : genes)
		{
			s += "\t" + gene;
		}
		return s;
	}
}
