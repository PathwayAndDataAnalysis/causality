package org.cbio.causality;

import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math.stat.regression.GLSMultipleLinearRegression;
import org.apache.commons.math.stat.regression.MultipleLinearRegression;
import org.cbio.causality.idmapping.EntrezGene;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class TempTest
{


	@Test
	@Ignore
	public void doMath()
	{
		GLSMultipleLinearRegression reg = new GLSMultipleLinearRegression();
		double[] Y = {3, 6, 8};
		double[][] X = {new double[]{3, 6, 7}, new double[]{6, 4, 7}, new double[]{6, 2, 9}};
		Covariance cov = new Covariance(X);
		reg.newSampleData(Y, X, cov.getCovarianceMatrix().getData());
		System.out.println(reg.estimateRegressionParameters());
	}

	@Test
	@Ignore
	public void findGeneSymbols()
	{
		String s = "Identifier\tDatabase\n" +
			"HMDB00122\tHMDB\n" +
			"ENSG00000113161\tEnsembl\n" +
			"ENSG00000129965\tEnsembl\n" +
			"ENSG00000066044\tEnsembl\n" +
			"ENSG00000141510\tEnsembl\n" +
			"ENSG00000133101\tEnsembl\n" +
			"ENSG00000145386\tEnsembl\n" +
			"ENSG00000134057\tEnsembl\n" +
			"1026\tEntrez Gene\n" +
			"ENSG00000120907\tEnsembl\n" +
			"ENSG00000170214\tEnsembl\n" +
			"WP534\tWikiPathways\n" +
			"WP534\tWikiPathways\n" +
			"WP143\tWikiPathways\n" +
			"ENSG00000105221\tEnsembl\n" +
			"ENSG00000117461\tEnsembl\n" +
			"ENSG00000110931\tEnsembl\n" +
			"ENSG00000006831\tEnsembl\n" +
			"ENSG00000108443\tEnsembl\n" +
			"1978\tEntrez Gene\n" +
			"53632\tEntrez Gene\n" +
			"ENSG00000104812\tEnsembl\n" +
			"ENSG00000103197\tEnsembl\n" +
			"5562\tEntrez Gene\n" +
			"ENSG00000142208\tEnsembl\n" +
			"31\tEntrez Gene\n" +
			"ENSG00000174697\tEnsembl\n" +
			"ENSG00000072310\tEnsembl\n" +
			"HMDB00058\tHMDB\n" +
			"ENSG00000171105\tEnsembl\n" +
			"ENSG00000078142\tEnsembl\n" +
			"51422\tEntrez Gene\n" +
			"2475\tEntrez Gene\n" +
			"5571\tEntrez Gene\n" +
			"ENSG00000145675\tEnsembl\n" +
			"ENSG00000125695\tEnsembl\n" +
			"5565\tEntrez Gene\n" +
			"ENSG00000165699\tEnsembl\n" +
			"ENSG00000051382\tEnsembl\n" +
			"ENSG00000121879\tEnsembl\n" +
			"ENSG00000082146\tEnsembl\n" +
			"ENSG00000159346\tEnsembl\n" +
			"ENSG00000155846\tEnsembl\n" +
			"ENSG00000079435\tEnsembl\n" +
			"5564\tEntrez Gene\n" +
			"ENSG00000182621\tEnsembl\n" +
			"ENSG00000064489\tEnsembl\n" +
			"ENSG00000103319\tEnsembl\n" +
			"1375\tEntrez Gene\n" +
			"6517\tEntrez Gene\n" +
			"ENSG00000169710\tEnsembl\n" +
			"ENSG00000175634\tEnsembl\n" +
			"1374\tEntrez Gene\n" +
			"200186\tEntrez Gene\n" +
			"HMDB01532\tHMDB\n" +
			"HMDB00464\tHMDB\n" +
			"ENSG00000101076\tEnsembl\n" +
			"ENSG00000116678\tEnsembl\n" +
			"5209\tEntrez Gene\n" +
			"ENSG00000105851\tEnsembl\n" +
			"HMDB03540\tHMDB\n" +
			"ENSG00000125520\tEnsembl\n" +
			"ENSG00000111713\tEnsembl\n" +
			"ENSG00000004660\tEnsembl\n" +
			"ENSG00000118046\tEnsembl\n" +
			"HMDB01175\tHMDB\n" +
			"1938\tEntrez Gene\n" +
			"51719\tEntrez Gene\n" +
			"126129\tEntrez Gene\n" +
			"5567\tEntrez Gene\n" +
			"ENSG00000171608\tEnsembl\n" +
			"5563\tEntrez Gene\n" +
			"57521\tEntrez Gene\n" +
			"ENSG00000105647\tEnsembl\n" +
			"ENSG00000181092\tEnsembl\n" +
			"5568\tEntrez Gene\n" +
			"32\tEntrez Gene";

		List<String> list = new ArrayList<String>();
		for (String s1 : s.split("\n"))
		{
			String[] tok = s1.split("\t");
			if (tok[1].equals("Entrez Gene"))
			{
				String symbol = EntrezGene.getSymbol(tok[0]);
				if (symbol != null) list.add(symbol);
			}
		}
		Collections.sort(list);
		for (String s1 : list)
		{
			System.out.print("\t" + s1);
		}
	}
}
