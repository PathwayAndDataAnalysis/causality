package org.cbio.causality.alteredreaction;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.AlterationProvider;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class AlteredReactionFinder
{
	public static void main(String[] args) throws IOException
	{
		AlterationProvider ap = new CBioPortalAccessor(PortalDataset.ENDOMETRIAL);
		AlteredReactionFinder arf = new AlteredReactionFinder();
		List<Reaction> list = arf.getRankedReactions(ap);
		for (int i = 0; i < 100; i++)
		{
			System.out.println(list.get(i));
		}
	}
	public List<Reaction> getRankedReactions(AlterationProvider ap)
	{
		List<Reaction> list = readReactions();
		fillAlterations(list, ap);
		removeSubsets(list);
		Collections.sort(list, new Comparator<Reaction>()
		{
			@Override
			public int compare(Reaction r1, Reaction r2)
			{
				return r2.coverage.compareTo(r1.coverage);
			}
		});
		return list;
	}

	public void removeSubsets(List<Reaction> list)
	{
		List<Reaction> remove = new ArrayList<Reaction>();

		for (Reaction r1 : list)
		{
			for (Reaction r2 : list)
			{
				if (r1 == r2) continue;

				if (r1.isSubsetOf(r2)) remove.add(r1);
			}
		}
		list.removeAll(remove);
	}

	public void fillAlterations(List<Reaction> list, AlterationProvider ap)
	{
		Set<String> notFound = new HashSet<String>();

		List<Reaction> remove = new ArrayList<Reaction>();
		for (Reaction r : list)
		{
			for (String gene : r.genes)
			{
				if (notFound.contains(gene)) continue;

				AlterationPack alt = ap.getAlterations(gene);
				if (alt == null)
				{
					notFound.add(gene);
				}
				else if (alt.isAltered())
				{
					r.alterations.put(gene, alt);
				}
			}
			r.genes.retainAll(r.alterations.keySet());
			if (r.genes.size() < 2) remove.add(r);
		}
		list.removeAll(remove);

		for (Reaction r : list)
		{
			r.sortGenes();
			r.fillChanges();
		}
	}

	public List<Reaction> readReactions()
	{
		List<Reaction> list = new ArrayList<Reaction>();

		Scanner scan = new Scanner(
			AlteredReactionFinder.class.getResourceAsStream("related-genes.txt"));

		while (scan.hasNextLine())
		{
			String line = scan.nextLine();

			String[] token = line.split("\t");

			Reaction r = new Reaction();
			r.ID = token[0];
			r.genes.addAll(Arrays.asList(token).subList(1, token.length));
			list.add(r);
		}
		scan.close();

		return list;
	}
}
