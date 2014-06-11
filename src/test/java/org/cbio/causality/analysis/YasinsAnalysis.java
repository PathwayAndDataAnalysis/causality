package org.cbio.causality.analysis;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.ProteinReference;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.Searcher;
import org.biopax.paxtools.pattern.constraint.ConBox;
import org.biopax.paxtools.pattern.constraint.Field;
import org.biopax.paxtools.pattern.constraint.OR;
import org.biopax.paxtools.pattern.constraint.PathConstraint;
import org.cbio.causality.util.Histogram;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class YasinsAnalysis
{
	public static void main(String[] args) throws IOException
	{
		YasinsAnalysis ya = new YasinsAnalysis();
		ya.run();
	}

	@Test
	@Ignore
	public void run() throws IOException
	{
		Set<String> select = getRPPAGenes();
		System.out.println("select.size() = " + select.size());

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Projects/biopax-pattern/All-Data.owl"));
		System.out.println("Model loaded");

		Map<String, Set<String>> gene2pubmed = new HashMap<String, Set<String>>();

		Pattern p = new Pattern(ProteinReference.class, "PR");
		p.add(ConBox.erToPE(), "PR", "SPE");
		p.add(ConBox.linkToComplex(), "SPE", "PE");
		p.add(ConBox.peToInter(), "PE", "Inter");
		p.add(new OR(new MappedConst(new PathConstraint("Interaction/xref:PublicationXref"), 0, 1),
			new MappedConst(new PathConstraint("Interaction/evidence/xref:PublicationXref"), 0, 1)),
			"Inter", "XR");
		p.add(new Field("Xref/db", Field.Operation.INTERSECT, "Pubmed"), "XR");

		List<Match> matches = Searcher.searchPlain(model, p);
		for (Match m : matches)
		{
			ProteinReference pr = (ProteinReference) m.get("PR", p);
			String sym = getSymbol(pr);
			if (sym != null)
			{
				Xref xr = (Xref) m.get("XR", p);
				if (!gene2pubmed.containsKey(sym)) gene2pubmed.put(sym, new HashSet<String>());
				gene2pubmed.get(sym).add(xr.getId());
			}
		}

		Histogram his0 = new Histogram(20);
		his0.setBorderAtZero(true);
		Histogram his1 = new Histogram(20);
		his1.setBorderAtZero(true);
		for (String gene : gene2pubmed.keySet())
		{
			Set<String> ids = gene2pubmed.get(gene);

			if (select.contains(gene)) his1.count(ids.size());
			else his0.count(ids.size());
		}
		his0.printDensity();
		his1.printDensity();

		int absent = 0;
		for (String gene : select)
		{
			if (!gene2pubmed.containsKey(gene)) absent++;
		}
		System.out.println("absent = " + absent);

		for (int i = 0; i < 10; i++)
		{
			int a = 2;
			System.out.println("a = " + a);
		}

		List<String> list = new ArrayList<String>(select);
		Collections.sort(list);
		BufferedWriter writer = new BufferedWriter(new FileWriter("gene_to_pubmed.txt"));
		writer.write("Symbol\t# of pubmed IDs\tpubmed IDs");
		for (String gene : list)
		{
			Set<String> ids = gene2pubmed.get(gene);
			writer.write("\n" + gene + "\t" + (ids == null ? "0" : ids.size()) + "\t");

			if (ids != null)
			{
				List<String> idlist = new ArrayList<String>(ids);
				Collections.sort(idlist);
				writer.write(idlist.toString());
			}
		}

		writer.close();
	}

	private String getSymbol(ProteinReference pr)
	{
		for (Xref xref : pr.getXref())
		{
			if (xref.getDb() != null && xref.getDb().equals("HGNC Symbol"))
			{
				return xref.getId();
			}
		}
		return null;
	}

	private Set<String> getRPPAGenes()
	{
		Set<String> set = new HashSet<String>();
		Scanner sc = new Scanner(YasinsAnalysis.class.getResourceAsStream("bp_pancan_input.txt"));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			Collections.addAll(set, line.split("\t")[2].split("\\|"));
		}

		return set;
	}

}
