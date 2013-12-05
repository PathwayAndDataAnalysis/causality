package org.cbio.causality.ubique;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.SmallMoleculeReference;
import org.biopax.paxtools.model.level3.XReferrable;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFSearcher;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.biopax.paxtools.pattern.util.HGNC;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

/**
 * @author Ozgun Babur
 */
public class SIFProducer
{
	public static void main(String[] args) throws FileNotFoundException
	{
		SimpleIOHandler reader = new SimpleIOHandler();
		Model model = reader.convertFromOWL(new FileInputStream(
			"/home/ozgun/Projects/biopax-pattern/All-Human-Data.owl"));

		SIFSearcher searcher = new SIFSearcher(new Fetcher(),
			SIFType.CONSUMPTION_CONTROLLED_BY, SIFType.CONTROLS_PRODUCTION_OF, SIFType.NEIGHBOR_OF);

		searcher.searchSIF(model, new FileOutputStream("sif-for-ubiques.txt"), false);
	}

	static class Fetcher implements IDFetcher
	{
		@Override
		public String fetchID(BioPAXElement ele)
		{
			if (ele instanceof SmallMoleculeReference)
			{
				return ChemicalNameNormalizer.nornmalize(
					((SmallMoleculeReference) ele).getDisplayName());
			}
			else if (ele instanceof XReferrable)
			{
				for (Xref xr : ((XReferrable) ele).getXref())
				{
					String db = xr.getDb();
					if (db != null)
					{
						db = db.toLowerCase();
						if (db.startsWith("hgnc"))
						{
							String id = xr.getId();
							if (id != null)
							{
								String symbol = HGNC.getSymbol(id);
								if (symbol != null && !symbol.isEmpty())
								{
									return symbol;
								}
							}
						}
					}
				}
			}

			return null;
		}
	}
}
