package org.cbio.causality.data.tcgafile;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import org.cbio.causality.util.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SimpleAnalysis
{
	private final static String TCGADIR = "/home/ozgun/Documents/TCGA/broad";

	//----- Expression plot -----------------

	public static void plotExpressions(String gene, String tcgaDir, String outFile) throws IOException
	{
		Set<String> genes = new HashSet<String>(Arrays.asList(gene));
		List<double[]> valList = new ArrayList<double[]>();
		List<String> names = new ArrayList<String>();

		File[] files = getCaseDirs(tcgaDir);

		for (File dir : files)
		{
			ExpressionReader er = new ExpressionReader(dir.getPath() + File.separator + "expression.txt", genes);
			Set<String> set = er.getSamples();
			double[] vals = er.getGeneAlterationArray(gene, set.toArray(new String[set.size()]));
			valList.add(vals);
			names.add(dir.getName());
		}

		ArrayUtil.prepareForBoxPlotR(valList, names, outFile);
	}

	public static void plotGeneExpression() throws IOException
	{
		String tcgaDir = "/home/ozgun/Documents/TCGA/";
		String gene = "NRAS";
		plotExpressions(gene, TCGADIR, tcgaDir + "/exp-dists/" + gene + ".txt");
	}

	//------- Expression correlation

	public static void reportCorrelations() throws FileNotFoundException, MathException
	{
		reportCorrelations("AGAP2", "CCND1", TCGADIR);
	}

	public static void reportCorrelations(String gene1, String gene2, String tcgaDir) throws FileNotFoundException, MathException
	{
		System.out.println(gene1 + " and " + gene2);
		Set<String> genes = new HashSet<String>(Arrays.asList(gene1, gene2));
		for (File dir : getCaseDirs(tcgaDir))
		{
			Object[] readers = prepareReaders(dir, genes);
			ExpressionReader er = (ExpressionReader) readers[0];
			CNAReader cr = (CNAReader) readers[1];
			MutationReader mr = (MutationReader) readers[2];

			String[] cases = getCommonCases(er, cr, mr);

			if (cases.length < 10) continue;

			Map<String, double[]> expMap = getUnalteredExpressions(genes, er, cr, mr, cases);

			double corr = Pearson.correlation(expMap.get(gene1), expMap.get(gene2));
			double pval = Pearson.corrPval(expMap.get(gene1), expMap.get(gene2));

			if (pval < 0.05)
			System.out.println(dir.getName() + " = " + corr + "\tpval = " + pval + "\tsss = " + cases.length);
		}
	}

	private static Map<String, double[]> getUnalteredExpressions(Set<String> genes,
		ExpressionReader er, CNAReader cr, MutationReader mr, String[] samples)
	{
		Map<String, double[]> map = new HashMap<String, double[]>();

		boolean[] b = getAltered(genes, cr, mr, samples);
		b = ArrayUtil.negate(b);

		for (String gene : genes)
		{
			double[] exp = er.getGeneAlterationArray(gene, samples);
			exp = ArrayUtil.subset(exp, b);
			map.put(gene, exp);
		}
		return map;
	}

	//--------- Alteration to expression ------------

	public static void plotEffectOfAlterationsOnExpression() throws FileNotFoundException
	{
		plotEffectOfAlterationsOnExpression("DLEU1", TCGADIR, "TFDP1", "CDK4", "RB1", "CDKN2A");
//		plotEffectOfAlterationsOnExpression(new String[]{"E2F1", "CDC14B", "CDC45L", "CDC6", "CDC7L1"},
//			TCGADIR + File.separator + "GBM",
//			"CDK4", "RB1", "CDKN2A");
	}

	public static void plotEffectOfAlterationsOnExpression(String target, String tcgaDir,
		String... alteredGenes) throws FileNotFoundException
	{
		for (File dir : getCaseDirs(tcgaDir))
		{
			reportExpressionChange(target, dir, alteredGenes);
		}
	}

	public static void plotEffectOfAlterationsOnExpression(String[] targets, String studyDir,
		String... alteredGenes) throws FileNotFoundException
	{
		for (String target : targets)
		{
			System.out.println("target = " + target);
			File dir = new File(studyDir);
			reportExpressionChange(target, dir, alteredGenes);
		}
	}

	private static void reportExpressionChange(String target, File dir, String[] alteredGenes) throws FileNotFoundException
	{
		Set<String> genes = new HashSet<String>(Arrays.asList(alteredGenes));
		genes.add(target);
		Object[] readers = prepareReaders(dir, genes);
		ExpressionReader er = (ExpressionReader) readers[0];
		CNAReader cr = (CNAReader) readers[1];
		MutationReader mr = (MutationReader) readers[2];

		String[] cases = getCommonCases(er, cr, mr);

		boolean[] b1 = getAltered(new HashSet<String>(Arrays.asList(alteredGenes)), cr, mr, cases);
//		boolean[] b1 = getAltered(new HashSet<String>(Arrays.asList(alteredGenes)), er, cr, mr, cases);
		boolean[] b2 = getAltered(Collections.singleton(target), cr, mr, cases);
//		boolean[] b2 = new boolean[b1.length];
		b2 = ArrayUtil.negate(b2);


		boolean[] test = ArrayUtil.getAND(b1, b2);
		boolean[] ctrl = ArrayUtil.getAND(ArrayUtil.negate(b1), b2);

		double[] exp = er.getGeneAlterationArray(target, cases);
		double[] testVals = ArrayUtil.subset(exp, test);
		double[] ctrlVals = ArrayUtil.subset(exp, ctrl);

		double pval = StudentsT.getPValOfMeanDifference(testVals, ctrlVals);
		if (pval < 0.05)
		{
			System.out.println(dir.getName() + "\tMean diff = " +
				(Summary.mean(testVals) - Summary.mean(ctrlVals)) + "\tpval = " + pval +
				"\tss= " + testVals.length + " " + ctrlVals.length);
		}
	}

	//--------- Compare two gene CNAs----------------------------------

	public static void compareTwoGeneCNAs() throws FileNotFoundException
	{
		String gene1 = "MDM2";
		String gene2 = "FRS2";

		final Map<String, Integer> cnt = new HashMap<String, Integer>();

		File dir = new File(TCGADIR);
		for (File studyDir : dir.listFiles())
		{
			CNAReader cr = new CNAReader(studyDir.getPath() + File.separator + "copynumber.txt",
				new HashSet<String>(Arrays.asList(gene1, gene2)), false, 2);
			Set<String> sampleSet = cr.getSamples();
			String[] samples = sampleSet.toArray(new String[sampleSet.size()]);
			int[] cna1 = cr.getGeneAlterationArray(gene1, samples);
			int[] cna2 = cr.getGeneAlterationArray(gene2, samples);

			for (int i = 0; i < samples.length; i++)
			{
				String s = (cna1[i] < 0 ? "" : " ") + cna1[i] + "\t" + (cna2[i] < 0 ? "" : " ") + cna2[i];

				if (cnt.containsKey(s)) cnt.put(s, cnt.get(s) + 1);
				else cnt.put(s, 1);
			}
		}

		List<String> keys = new ArrayList<String>(cnt.keySet());
		Collections.sort(keys, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return cnt.get(o2).compareTo(cnt.get(o1));
			}
		});

		for (String s : keys)
		{
			System.out.println(s + "\t" + cnt.get(s));
		}
	}


	//--------- Best correlating alteration to expression -------------

	public static void findBestCorrelatingAlterationsOfExpression(String target, String studyDir,
		double expThr, boolean lowExp) throws FileNotFoundException
	{
		Object[] r = prepareReaders(new File(studyDir), null);
		String[] cases = getCommonCases(r);


	}

	//--------- List mutation types ----------

	public static void listMutationTypes() throws FileNotFoundException
	{
		String gene = "NCOR2";
		String file = TCGADIR + File.separator + "SARC/mutation.maf";
		TermCounter tc = new TermCounter();
		TermCounter t2 = new TermCounter();

		System.out.println("gene = " + gene);
		Scanner sc = new Scanner(new File(file));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (!line.startsWith(gene)) continue;
			String id = line.substring(0, line.indexOf("\t"));
			if (!id.equals(gene)) continue;
			String[] token = line.split("\t");
			String type = token[8];
			String ch = token[49];
			if (!type.isEmpty()) tc.addTerm(type);
			if (!ch.isEmpty()) t2.addTerm(ch);
		}
		tc.print();
		t2.print();
	}

	//--------- List alteration frequencies ---------

	public static void listAlterationFrequencies() throws FileNotFoundException
	{
		String gene = "TP53";
		String base = TCGADIR + File.separator + "SARC/";
		MutationReader mr = new MutationReader(base + "mutation.maf");
		CNAReader cr = new CNAReader(base + "copynumber.txt", Collections.singleton(gene));
		Set<String> cases = mr.getSamples();
		cases.retainAll(cr.getSamples());
		String[] samples = cases.toArray(new String[cases.size()]);
		boolean[] mut = mr.getGeneAlterationArray(gene, samples);
		int[] cna = cr.getGeneAlterationArray(gene, samples);

		int mCnt = 0;
		int ampCnt = 0;
		int delCnt = 0;
		int mutOrDel = 0;
		for (int i = 0; i < samples.length; i++)
		{
			if (mut[i]) mCnt++;
			if (cna[i] > 0) ampCnt++;
			else if (cna[i] < 0) delCnt++;

			if (mut[i] || cna[i] < 0) mutOrDel++;
		}
		System.out.println("samples = " + samples.length);
		System.out.println("mCnt = " + mCnt);
		System.out.println("ampCnt = " + ampCnt);
		System.out.println("delCnt = " + delCnt);
		System.out.println("mutOrDel = " + mutOrDel);
	}

	//--------- Commons -----------

	private static File[] getCaseDirs(String tcgaDir)
	{
		File[] files = new File(tcgaDir).listFiles();
		Arrays.sort(files, new Comparator<File>()
		{
			@Override
			public int compare(File o1, File o2)
			{
				return o1.getName().compareTo(o2.getName());
			}
		});
		return files;
	}

	private static boolean[] getAltered(Set<String> genes, CNAReader cr, MutationReader mr, String[] samples)
	{
		boolean[] b = new boolean[samples.length];
		Arrays.fill(b, false);

		for (String gene : genes)
		{
			int[] cna = cr.getGeneAlterationArray(gene, samples);
			for (int i = 0; i < cna.length; i++)
			{
				if (cna[i] != 0) b[i] = true;
			}
			boolean[] mut = mr.getGeneAlterationArray(gene, samples);
			for (int i = 0; i < mut.length; i++)
			{
				if (mut[i]) b[i] = true;
			}
		}
		return b;
	}

	private static boolean[] getAltered(Set<String> genes, ExpressionReader er, CNAReader cr, MutationReader mr, String[] samples)
	{
		boolean[] b = new boolean[samples.length];
		Arrays.fill(b, false);

		for (String gene : genes)
		{
			int[] cna = cr.getExpVerifiedCNA(gene, samples, er.getGeneAlterationArray(gene, samples), 0.05);
			if (cna != null)
			{
				for (int i = 0; i < cna.length; i++)
				{
					if (cna[i] != 0) b[i] = true;
				}
			}
			boolean[] mut = mr.getGeneAlterationArray(gene, samples);
			for (int i = 0; i < mut.length; i++)
			{
				if (mut[i]) b[i] = true;
			}
		}
		return b;
	}

	private static Object[] prepareReaders(File dir, Set<String> genes) throws FileNotFoundException
	{
		ExpressionReader er = new ExpressionReader(dir.getPath() + File.separator + "expression.txt", genes);
		CNAReader cr = new CNAReader(dir.getPath() + File.separator + "copynumber.txt", genes);
		MutationReader mr = new MutationReader(dir.getPath() + File.separator + "mutation.maf");
		return new Object[]{er, cr, mr};
	}

	private static String[] getCommonCases(Object[] r)
	{
		return getCommonCases((ExpressionReader) r[0], (CNAReader) r[1], (MutationReader) r[2]);
	}

	private static String[] getCommonCases(ExpressionReader er, CNAReader cr, MutationReader mr)
	{
		Set<String> erSet = er.getSamples();
		Set<String> crSet = cr.getSamples();
		Set<String> mrSet = mr.getSamples();

		Set<String> sampleSet = new HashSet<String>(erSet);
		sampleSet.retainAll(crSet);
		sampleSet.retainAll(mrSet);
		return sampleSet.toArray(new String[sampleSet.size()]);
	}

	//------ Main ---------------

	public static void main(String[] args) throws IOException, MathException
	{
//		plotGeneExpression();
//		reportCorrelations();
//		plotEffectOfAlterationsOnExpression();
//		listMutationTypes();
//		listAlterationFrequencies();
		compareTwoGeneCNAs();
	}
}
