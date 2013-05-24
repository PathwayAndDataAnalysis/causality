package org.cbio.causality.alteredreaction;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Interaction;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.AlterationProvider;
import org.cbio.causality.util.ModelExciser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class AlteredReactionFinder
{
	private static Log log = LogFactory.getLog(AlteredReactionFinder.class);

	public static void main(String[] args) throws IOException
	{
		AlterationProvider ap = new CBioPortalAccessor(PortalDataset.ENDOMETRIAL_MUT_CN);
		AlteredReactionFinder arf = new AlteredReactionFinder();
		String inFile = "/home/ozgun/Desktop/PC.owl";
		String outFile = "/home/ozgun/Desktop/reactions.owl";
		arf.exciseAlteredReactions(ap, 500, inFile, outFile);
	}

	public void exciseAlteredReactions(AlterationProvider ap, int maxSizeOfResults, String inFile,
		String outFile)
	{
		SimpleIOHandler handler = new SimpleIOHandler();
		try
		{
			List<Reaction> list = getRankedReactions(ap);

			Model model = handler.convertFromOWL(new FileInputStream(inFile));

			List<ModelExciser.PathwayTicket> tickets = new ArrayList<ModelExciser.PathwayTicket>();

			int i = 0;
			for (Reaction r : list)
			{
				if (i++ == maxSizeOfResults) break;

				Interaction inter = (Interaction) model.getByID(r.ID);

				if (inter == null)
				{
					log.error("Cannot find reaction " + r.ID + " in the model " + inFile);
					continue;
				}

				tickets.add(new ModelExciser.PathwayTicket(r.getGeneNames(),
					"Coverage: " + r.coverage, inter));
			}

			log.info("Creating " + tickets.size() + " pathways.");

			model = ModelExciser.excise(model, tickets, true, true);

			handler.convertToOWL(model, new FileOutputStream(outFile));
		}
		catch (FileNotFoundException e)
		{
			log.error("Cannot excise.", e);
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
				return r2.score.compareTo(r1.score);
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

				if (r1.isSubsetOf(r2) && r1.score <= r2.score)
				{
					remove.add(r1);
				}
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
