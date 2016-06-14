package org.cbio.causality.analysis;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.PatternBox;
import org.biopax.paxtools.pattern.Searcher;
import org.biopax.paxtools.pattern.constraint.PEChainsIntersect;
import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;
import org.cbio.causality.util.TermCounter;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class AlexAnalysis
{
	public static void main(String[] args) throws FileNotFoundException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream(
			"/home/ozgun/Projects/biopax-pattern/All-Data.owl"));

		TermCounter tc = new TermCounter();
		for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class))
		{
			for (EntityFeature f : pe.getFeature())
			{
				if (f instanceof ModificationFeature)
				{
					SequenceModificationVocabulary type =
						((ModificationFeature) f).getModificationType();

					if (type != null)
					{
						for (String term : type.getTerm())
						{
							if (term.contains("hosp")) tc.addTerm(term);
						}
					}
				}
			}
		}

		tc.print();
	}

	@Test
	@Ignore
	public void generatePhosphoGraph() throws IOException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream(
			"/home/ozgun/Projects/biopax-pattern/All-Data.owl"));

		Pattern p = PatternBox.controlsStateChange();
		p.add(new PEChainsIntersect(false),
			"input simple PE", "input PE", "output simple PE", "output PE");
		List<Match> matches = Searcher.searchPlain(model, p);

		Map<String, Set<String>> mediators = new HashMap<String, Set<String>>();
		Map<String, Set<String>> gainedFeats = new HashMap<String, Set<String>>();
		Map<String, Set<String>> lostFeats = new HashMap<String, Set<String>>();

		for (Match match : matches)
		{
			Protein pIn = (Protein) match.get("input simple PE", p);
			Protein pOut = (Protein) match.get("output simple PE", p);
			ProteinReference changedPr = (ProteinReference) match.get("changed PR", p);
			ProteinReference ctrlPR = (ProteinReference) match.get("controller PR", p);
			Protein controller = (Protein) match.get("controller simple PE", p);
			Control ctrl = (Control) match.get("Control", p);
			Conversion cnv = (Conversion) match.get("Conversion", p);

			String source = getSymbol(ctrlPR);
			String target = getSymbol(changedPr);

			if (source == null || target == null) continue;

			Set<ModificationFeature>[] mods =
				DifferentialModificationUtil.getChangedModifications(pIn, pOut);

			if (mods[0].isEmpty() && mods[1].isEmpty()) continue;

			Set<String> gained = getPhosphoStrings(mods[0], changedPr);
			Set<String> lost = getPhosphoStrings(mods[1], changedPr);

			if (gained.isEmpty() && lost.isEmpty()) continue;

			boolean inactive = hasInactiveTag(controller);
			boolean inhibiting = ctrl.getControlType() != null &&
				ctrl.getControlType().name().startsWith("I");

			if ((inactive || inhibiting) && !(inactive && inhibiting)) // xor
			{
				Set<String> temp = gained;
				gained = lost;
				lost = temp;
			}

			String key = source + "\t" + SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag() + "\t" + target;

			if (!mediators.containsKey(key)) mediators.put(key, new HashSet<String>());
			if (!gainedFeats.containsKey(key)) gainedFeats.put(key, new HashSet<String>());
			if (!lostFeats.containsKey(key)) lostFeats.put(key, new HashSet<String>());

			mediators.get(key).add(cnv.getUri());
			mediators.get(key).add(ctrl.getUri());

			gainedFeats.get(key).addAll(gained);
			lostFeats.get(key).addAll(lost);
		}

		BufferedWriter writer = new BufferedWriter(new FileWriter(
			"/home/ozgun/Desktop/phospho.sif"));

		for (String key : mediators.keySet())
		{
			writer.write(key + "\t");

			writer.write(concat(mediators.get(key)) + "\t");
			writer.write(concat(sortModif(gainedFeats.get(key))) + "\t");
			writer.write(concat(sortModif(lostFeats.get(key))) + "\n");
		}


		writer.close();
	}

	private List<String> sortModif(Set<String> modifs)
	{
		List<String> list = new ArrayList<String>(modifs);
		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				if (o1.length() < 2) return -1;
				if (o2.length() < 2) return 1;

				return new Integer(o1.substring(1)).compareTo(new Integer(o2.substring(1)));
			}
		});

		return list;
	}

	private String concat(Collection<String> ss)
	{
		String x = "";
		boolean first = true;
		for (String s : ss)
		{
			if (first) first = false;
			else x += " ";

			x += s;
		}
		return x;
	}

	private void writeModif(Writer writer, List<String> list) throws IOException
	{
		writer.write("\t");

		for (int i = 0; i < list.size(); i++)
		{
			if (i > 0) writer.write(" ");
			writer.write(list.get(i));
		}
	}

	private boolean hasInactiveTag(Protein p)
	{
		for (EntityFeature f : p.getFeature())
		{
			if (f instanceof ModificationFeature)
			{
				SequenceModificationVocabulary type =
					((ModificationFeature) f).getModificationType();

				if (type != null)
				{
					for (String term : type.getTerm())
					{
						if (term.contains("inactive")) return true;
					}
				}
			}
		}
		return false;
	}

	private String getSymbol(EntityReference er)
	{
		for (Xref xref : er.getXref())
		{
			if (xref.getDb() != null && xref.getDb().equals("HGNC Symbol")) return xref.getId();
		}
		return null;
	}

	private Set<String> getPhosphoStrings(Collection<ModificationFeature> mfs, ProteinReference pr)
	{
		Set<String> set = new HashSet<String>();
		for (ModificationFeature mf : mfs)
		{
			String ps = getPhoshoString(mf, pr);
			if (ps != null) set.add(ps);
		}

		return set;
	}

	private String getPhoshoString(ModificationFeature mf, ProteinReference pr)
	{
		if (!isPhosphorylation(mf)) return null;

		String s;

		if (contains(mf, "serine"))
		{
			s = "S";
		}
		else if (contains(mf, "threonine"))
		{
			s = "T";
		}
		else if (contains(mf, "tyrosine"))
		{
			s = "Y";
		}
		else if (contains(mf, "lysine"))
		{
			s = "K";
		}
		else if (contains(mf, "cysteine"))
		{
			s = "C";
		}
		else
		{
			s = "P";
		}

		SequenceLocation loc = mf.getFeatureLocation();
		if (loc instanceof SequenceSite)
		{
			int site = ((SequenceSite) loc).getSequencePosition();

			s += site;
		}

		return s;
	}

	private boolean contains(ModificationFeature mf, String s)
	{
		SequenceModificationVocabulary type = mf.getModificationType();
		if (type == null) return false;

		for (String term : type.getTerm())
		{
			if (term.toLowerCase().contains(s)) return true;
		}
		return false;
	}

	private static final Set<String> TERMS = new HashSet<String>(Arrays.asList(
		"O-phospho-L-serine",
		"O4'-phospho-L-tyrosine",
		"O-phospho-L-threonine",
		"phosphorylated residue",
		"N6-pyridoxal phosphate-L-lysine",
		"phosphorylation",
		"O-phosphopantetheine-L-serine",
		"S-phospho-L-cysteine"));

	private boolean isPhosphorylation(ModificationFeature mf)
	{
		SequenceModificationVocabulary type = mf.getModificationType();
		if (type != null)
		{
			for (String term : type.getTerm())
			{
				if (TERMS.contains(term)) return true;
			}
		}
		return false;
	}
}
