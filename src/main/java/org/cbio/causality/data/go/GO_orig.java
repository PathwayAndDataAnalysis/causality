package org.cbio.causality.data.go;

import org.cbio.causality.util.FishersExactTest;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * @author Ozgun Babur
 */
public class GO_orig
{
	private static Map<String, Term> terms = new HashMap<String, Term>();
	private static Map<String, Set<Term>> gene2Term = new HashMap<String, Set<Term>>();

	public static Map<String, Double> getEnrichedTerms(Set<String> selectedGenes,
		Set<String> backgroundGenes)
	{
		if (!backgroundGenes.containsAll(selectedGenes)) throw new IllegalArgumentException(
			"Background genes have to contain all the selected genes.");

		Map<Term, Integer> selectionCnt = count(selectedGenes);
		Map<Term, Integer> backgroundCnt = count(backgroundGenes);

//		Map<Integer, Integer> bgCounts = getBgCounts(backgroundCnt);

		Map<String, Double> map = new HashMap<String, Double>();

		for (Term term : selectionCnt.keySet())
		{
			double pval = FishersExactTest.calcEnrichmentPval(backgroundGenes.size(),
				backgroundCnt.get(term), selectedGenes.size(), selectionCnt.get(term));

//			if (bgCounts.containsKey(selectionCnt.get(term)))
//			{
//				pval *= bgCounts.get(selectionCnt.get(term));
//			}

			map.put(term.name, pval);
		}

		return map;
	}

	private static Map<Integer, Integer> getBgCounts(Map<Term, Integer> map)
	{
		Map<Integer, Integer> cnt = new HashMap<Integer, Integer>();

		for (Term term : map.keySet())
		{
			Integer i = map.get(term);

			for (int j = 1; j <= i; j++)
			{
				if (!cnt.containsKey(j)) cnt.put(j, 1);
				else cnt.put(j, cnt.get(j) + 1);
			}
		}

		return cnt;
	}

	private static Map<Term, Integer> count(Set<String> genes)
	{
		Map<Term, Integer> cnt = new HashMap<Term, Integer>();

		for (String gene : genes)
		{
			if (gene2Term.containsKey(gene))
			{
				for (Term term : gene2Term.get(gene))
				{
					if (!cnt.containsKey(term)) cnt.put(term, 1);
					else cnt.put(term, cnt.get(term) + 1);
				}
			}
		}

		return cnt;
	}

	private static void readHumanAnnotation() throws IOException
	{

		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
			GO_orig.class.getResourceAsStream("gene_association.goa_human.gz"))));

		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			if (line.startsWith("!")) continue;

			String[] token = line.split("\t");

			if (token.length > 4)
			{
				if (!gene2Term.containsKey(token[2])) gene2Term.put(token[2], new HashSet<Term>());
				assert terms.containsKey(token[4]);
				gene2Term.get(token[2]).add(terms.get(token[4]));
			}
		}
	}

	private static void readOntology() throws IOException
	{
		Map<String, Set<String>> parents = new HashMap<String, Set<String>>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
			GO_orig.class.getResourceAsStream("gene_ontology.1_2.obo.txt.gz"))));

		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			if (line.startsWith("[Term]"))
			{
				Term term = new Term();

				line = reader.readLine();
				while (! (line == null || line.isEmpty()))
				{
					if (line.startsWith("id:")) term.id = getValue(line);
					else if (line.startsWith("name:")) term.name = getValue(line);
					else if (line.startsWith("namespace:"))
						term.namespace = Namespace.valueOf(getValue(line));
					else if (line.startsWith("is_a:"))
					{
						String id = line.split(" ")[1];
						if (!parents.containsKey(term.id)) parents.put(term.id, new HashSet<String>());
						parents.get(term.id).add(id);
					}

					line = reader.readLine();
				}

				terms.put(term.id, term);
			}
		}

		reader.close();


		for (String id : parents.keySet())
		{
			Term term = terms.get(id);

			for (String parent : parents.get(id))
			{
				assert terms.containsKey(parent);
				term.parents.add(terms.get(parent));
			}
		}
	}

	private static String getValue(String line)
	{
		return line.substring(line.indexOf(": ") + 2);
	}

	static
	{
		try
		{
			readOntology();
			readHumanAnnotation();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	static class Term
	{
		String id;
		String name;
		Namespace namespace;
		Set<Term> parents = new HashSet<Term>();

		@Override
		public int hashCode()
		{
			return id.hashCode();
		}

		@Override
		public boolean equals(Object obj)
		{
			return obj instanceof Term && id.equals(((Term) obj).id);
		}
	}

	enum Namespace
	{
		molecular_function,
		biological_process,
		cellular_component
	}

	public static void main(String[] args)
	{
		GO_orig go = new GO_orig();
	}
}
